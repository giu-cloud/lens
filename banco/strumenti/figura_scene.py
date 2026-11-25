"""Scene della cassaforte: radar (vento vero) | ECMWF letto dal mare | ECMWF corretto con gli analoghi. Nodi.

  VERO=<radar_vero.npy> python3 figura_scene.py <tabella.npz> <previsioni.npz> <cartella del metro> <uscita.png> area:n ...
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confronto_metodi import AREE, KN, carica


def main():
    tab, prev, metro, usc = sys.argv[1:5]
    richieste = [a.split(":") for a in sys.argv[5:]]
    T = carica(tab)
    P = np.load(prev)
    righe = []
    for area, n in richieste:
        ia = AREE.index(area)
        m = (T["area"] == ia) & T["cassa"] & np.isfinite(P["analoghi"])
        sc, cnt = np.unique(T["scena"][m], return_counts=True)
        S0 = np.array([T["S0"][m][T["scena"][m] == s][0] for s in sc])
        buone = sc[(cnt > np.percentile(cnt, 50)) & (S0 > 6)]
        rng = np.random.default_rng(1)
        for s in rng.choice(buone, min(int(n), len(buone)), replace=False):
            righe.append((area, s))
    fig, axs = plt.subplots(len(righe), 3, figsize=(10, 2.6 * len(righe)), layout="constrained")
    for r, (area, s) in enumerate(righe):
        z = np.load(f"{metro}/campi_{area}.npz")
        ny, nx = len(z["lat"]), len(z["lon"])
        m = T["scena"] == s
        j = T["j"][m].astype(int); i = T["i"][m].astype(int)
        mapp = []
        for val in (T["radar"][m], T["mod"][m], T["mod"][m] * P["analoghi"][m]):
            G = np.full((ny, nx), np.nan); G[j, i] = val * KN; mapp.append(G)
        vmax = np.nanpercentile(mapp[0], 98)
        tt = np.datetime64(int(T["t"][m][0]), "s")
        for c, (G, tit) in enumerate(zip(mapp, ["radar (vento vero)", "ECMWF letto dal mare", "ECMWF + analoghi"])):
            ax = axs[r, c]
            im = ax.pcolormesh(z["lon"], z["lat"], G, cmap="viridis", vmin=0, vmax=vmax, shading="auto")
            ax.contour(z["lon"], z["lat"], z["dist"], [0.01], colors="w", linewidths=0.6)
            ax.set_aspect(1 / np.cos(np.radians(z["lat"].mean()))); ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(tit, fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{area}\n{str(tt)[:16].replace('T', ' ')} UTC", fontsize=9)
        plt.colorbar(im, ax=axs[r, :], shrink=0.95, label="nodi", pad=0.01)
    plt.savefig(usc, dpi=75)


if __name__ == "__main__":
    main()
