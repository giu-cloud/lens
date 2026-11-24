"""Carte dell'atlante imparato dal radar (vento vero ÷ modello letto dal mare) per area e ottante di provenienza.

  python3 carte_atlante.py <tabella.npz> <radar_vero.npy> <cartella del metro> <area> <uscita.png>
"""
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AREE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]
OTT = ["N", "NE", "E", "SE", "S", "SO", "O", "NO"]


def main():
    tab, vero, metro, area, usc = sys.argv[1:6]
    T = np.load(tab)
    ia = AREE.index(area)
    m = T["area"] == ia
    R = np.load(vero)[m]; M = T["mod"][m]; j = T["j"][m].astype(int); i = T["i"][m].astype(int)
    ott = (np.round(T["dir"][m] / 45.0).astype(int) % 8)
    sc = T["scena"][m]
    z = np.load(f"{metro}/campi_{area}.npz")
    ny, nx = len(z["lat"]), len(z["lon"])
    fig, axs = plt.subplots(2, 4, figsize=(17, 8.2))
    for o in range(8):
        s = (ott == o) & (M > 4)
        k = j[s] * nx + i[s]
        num = np.bincount(k, weights=(R * M)[s], minlength=ny * nx); den = np.bincount(k, weights=(M * M)[s], minlength=ny * nx)
        n = np.bincount(k, minlength=ny * nx)
        C = np.where(n >= 8, num / np.maximum(den, 1e-9), np.nan).reshape(ny, nx)
        ax = axs.flat[o]
        im = ax.pcolormesh(z["lon"], z["lat"], C, cmap="RdBu_r", vmin=0.7, vmax=1.3, shading="auto")
        ax.contour(z["lon"], z["lat"], z["dist"], [0.01], colors="k", linewidths=0.6)
        ax.set_aspect(1 / np.cos(np.radians(z["lat"].mean())))
        ax.set_title(f"da {OTT[o]}: {len(np.unique(sc[s]))} passaggi oltre 8 nodi", fontsize=9)
        ax.tick_params(labelsize=7)
    fig.colorbar(im, ax=axs, shrink=0.7, label="vento vero (radar) ÷ modello letto dal mare")
    fig.suptitle(f"{area}: che cosa insegna il radar, per direzione di provenienza (blu = il modello ne ha troppo, rosso = troppo poco)", fontsize=11)
    plt.savefig(usc, dpi=70)


if __name__ == "__main__":
    main()
