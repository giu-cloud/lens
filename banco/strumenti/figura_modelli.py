"""Scene della cassaforte: radar (vento vero) | ECMWF | ICON-2I | insieme corretto (analoghi su tutti e due). Nodi.

  VERO=<radar_vero_icon.npy> python3 figura_modelli.py <tabella_ecmwf_sub> <tabella_icon> <prev_ecmwf> <prev_icon> <metro> <uscita.png> area:n ...
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
    te, ti, pe, pi, metro, usc = sys.argv[1:7]
    E = carica(te); I = np.load(ti)
    PE = np.load(pe)["analoghi"]; PI = np.load(pi)["analoghi"]
    righe = []
    for a in sys.argv[7:]:
        area, n = a.split(":")
        ia = AREE.index(area)
        m = (E["area"] == ia) & E["cassa"] & np.isfinite(PE) & np.isfinite(PI)
        sc, cnt = np.unique(E["scena"][m], return_counts=True)
        S0 = np.array([E["S0"][m][E["scena"][m] == s][0] for s in sc])
        buone = sc[(cnt > np.percentile(cnt, 50)) & (S0 > 6)]
        for s in np.random.default_rng(3).choice(buone, min(int(n), len(buone)), replace=False):
            righe.append((area, s))
    tit = ["radar (vento vero)", "ECMWF letto dal mare", "ICON-2I", "insieme corretto"]
    fig, axs = plt.subplots(len(righe), 4, figsize=(13, 2.8 * len(righe)), layout="constrained")
    axs = np.atleast_2d(axs)
    for r, (area, s) in enumerate(righe):
        z = np.load(f"{metro}/campi_{area}.npz"); ny, nx = len(z["lat"]), len(z["lon"])
        m = E["scena"] == s
        j = E["j"][m].astype(int); i = E["i"][m].astype(int)
        vals = [E["radar"][m], E["mod"][m], I["mod"][m], 0.65 * E["mod"][m] * PE[m] + 0.35 * I["mod"][m] * PI[m]]
        vmax = np.nanpercentile(vals[0] * KN, 98)
        tt = np.datetime64(int(E["t"][m][0]), "s")
        for c, v in enumerate(vals):
            G = np.full((ny, nx), np.nan); G[j, i] = v * KN
            ax = axs[r, c]
            im = ax.pcolormesh(z["lon"], z["lat"], G, cmap="viridis", vmin=0, vmax=vmax, shading="auto")
            ax.contour(z["lon"], z["lat"], z["dist"], [0.01], colors="w", linewidths=0.6)
            ax.set_aspect(1 / np.cos(np.radians(z["lat"].mean()))); ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(tit[c], fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{area}\n{str(tt)[:16].replace('T', ' ')} UTC", fontsize=9)
        plt.colorbar(im, ax=axs[r, :], shrink=0.95, label="nodi", pad=0.01)
    plt.savefig(usc, dpi=72)


if __name__ == "__main__":
    main()
