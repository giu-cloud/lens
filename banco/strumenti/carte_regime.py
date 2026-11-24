"""Carte medie del radar contro il modello, per tipo di vento, su un'area: dove il modello sottovaluta o sopravvaluta.

  python3 carte_regime.py <cartella dei costa_*.npz> <rilievo .rgr> <uscita> area [cartella del metro]

Senza la cartella del metro usa il vento ECMWF che accompagna ogni passaggio nel prodotto del radar (tutte le scene,
lettura del prodotto). Con la cartella del metro (campi_*.npz e sho_*.bin di costa_ecmwf.py) usa ECMWF letto come
Shift 0.54, dai soli punti di mare: è il confronto che conta per Lens. Per ogni scena la
direzione e la forza si prendono al largo (mediana delle celle oltre 8 km dalla terra); si tengono le scene con oltre
8 nodi al largo. Per ogni gruppo di direzioni: rapporto mediano radar ÷ modello cella per cella, solo dove ci sono
almeno 5 scene; sotto, il rilievo di Shift.
"""
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "costa"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_analisi as ca
from rilievo import Rilievo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

KN = 1.943844
GRUPPI = [("ponente e maestrale (da O a NO)", 270, 315), ("libeccio (da SO a OSO)", 225, 260),
          ("scirocco (da SE a S)", 135, 180), ("levante (da ENE a ESE)", 67.5, 112.5),
          ("tramontana e grecale (da NNO a NE)", 337.5, 45)]


def dentro(d, a, b):
    return (d >= a) & (d <= b) if a <= b else (d >= a) | (d <= b)


def main():
    dati, ril_p, usc, area = sys.argv[1:5]
    metro = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(usc, exist_ok=True)
    if metro:
        from diagnosi_rilievo import campi
        z, mu, mv = campi(metro, area)
        lat, lon = z["lat"], z["lon"]
        R = z["radar"].astype("f8"); M = np.hypot(mu, mv)
        lettura = "ECMWF 0,25 letto dal mare (Shift 0.54)"
    else:
        lat, lon, scene, _ = ca.carica(dati, area)
        R = np.array([np.hypot(s["su"], s["sv"]) for s in scene]); M = np.array([np.hypot(s["mu"], s["mv"]) for s in scene])
        mu = np.array([s["mu"] for s in scene]); mv = np.array([s["mv"] for s in scene])
        lettura = "ECMWF del prodotto del radar"
    mare = np.isfinite(R).any(axis=0)
    kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
    D = ndimage.distance_transform_edt(mare, sampling=(kmy, kmx))
    larg = (D > 8)[None] & np.isfinite(M)
    U = np.array([np.nanmedian(M[i][larg[i]]) if larg[i].any() else np.nan for i in range(len(R))])
    dd = np.array([(np.degrees(np.arctan2(-np.nanmedian(mu[i][larg[i]]), -np.nanmedian(mv[i][larg[i]]))) + 360) % 360
                   if larg[i].any() else np.nan for i in range(len(R))])
    rl, ro, H = Rilievo(ril_p).griglia(0, lat.min() - 0.02, lat.max() + 0.02, lon.min() - 0.02, lon.max() + 0.02)
    gruppi = [(g, a, b, (U * KN > 8) & dentro(np.nan_to_num(dd, nan=-999), a, b)) for g, a, b in GRUPPI]
    gruppi = [x for x in gruppi if x[3].sum() >= 6]
    if not gruppi:
        print(f"{area}: nessun tipo di vento con almeno 6 passaggi oltre 8 nodi"); return
    fig, axs = plt.subplots(1, len(gruppi), figsize=(4.6 * len(gruppi), 4.4), squeeze=False)
    norm = TwoSlopeNorm(vmin=0.6, vcenter=1.0, vmax=1.6)
    for ax, (g, a, b, sel) in zip(axs[0], gruppi):
        r = R[sel] / np.where(M[sel] > 1.5, M[sel], np.nan)
        n = np.isfinite(r).sum(axis=0)
        med = np.where(n >= 5, np.nanmedian(r, axis=0), np.nan)
        ax.contourf(ro, rl, np.where(H > 0, H, np.nan), levels=[0, 100, 300, 600, 900, 1300], cmap="Greys", alpha=0.55)
        ax.contour(ro, rl, H, levels=[1], colors="k", linewidths=0.6)
        im = ax.pcolormesh(lon, lat, med, cmap="RdBu_r", norm=norm, shading="nearest")
        ax.set_title(f"{g}\n{sel.sum()} passaggi", fontsize=9)
        ax.set_aspect(1 / np.cos(np.radians(lat.mean()))); ax.tick_params(labelsize=7)
    cb = fig.colorbar(im, ax=axs[0].tolist(), shrink=0.85, pad=0.02)
    cb.set_label("radar ÷ modello (rosso: il radar vede più vento)", fontsize=8)
    fig.suptitle(f"{area}: radar dei Sentinel-1 contro {lettura}, mediana cella per cella, vento oltre 8 nodi al largo", fontsize=10)
    dst = os.path.join(usc, f"regimi_{area}{'_mare' if metro else ''}.png")
    fig.savefig(dst, dpi=130, bbox_inches="tight")
    print(dst)


if __name__ == "__main__":
    main()
