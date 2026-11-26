"""Carte composite dello Stretto: vento ÷ modello, media sui passaggi del radar con il regime (N, S), oltre 4 m/s al largo.
  python3 figura_stretto.py <carte.npz> <atlante .lna> <png>"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
from controlla_atlante import leggi as leggi_lna
Z = np.load(sys.argv[1]); A = leggi_lna(sys.argv[2]); png = sys.argv[3]
S, N, O, E = A["dominio"]; h = A["h"]; NY, NX = h.shape
la_t = N - (np.arange(NY) + 0.5) * (N - S) / NY; lo_t = O + (np.arange(NX) + 0.5) * (E - O) / NX
lat, lon = Z["lat"], Z["lon"]
col = [("radar", "radar (vero), tutti i passaggi"), ("A", "A: primo Lens"), ("C", "C: A + livello del radar a 30 km"),
       ("B", "B: radar cella per cella"), ("CFD", "CFD neutra a 150 m")]
fig, ax = plt.subplots(2, 5, figsize=(26, 13))
norm = TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=2.0)
for r, (rg, nome_rg) in enumerate((("N", "vento da nord (tramontana, grecale)"), ("S", "vento da sud (scirocco, ostro)"))):
    for c, (k, tit) in enumerate(col):
        a = ax[r, c]
        im = a.pcolormesh(lon, lat, Z[f"{rg}_{k}"], cmap="RdBu_r", norm=norm, shading="nearest")
        a.contour(lo_t, la_t, (h > 0).astype(float), [0.5], colors="k", linewidths=0.6)
        a.contourf(lo_t, la_t, (h > 0).astype(float), [0.5, 1.5], colors=["#d9d4c7"])
        n = int(np.nanmax(Z[f"{rg}_{'radar' if k in ('radar', 'A', 'CFD') else 'B'}_n"]))
        a.set_title(f"{nome_rg}\n{tit}" + ("  (solo cassaforte)" if k in ("B", "C") else ""), fontsize=10)
        a.set_xlim(15.40, 15.75); a.set_ylim(37.90, 38.35); a.set_aspect(1 / np.cos(np.radians(38.1)))
fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.6, label="vento ÷ modello (media dei passaggi)")
plt.savefig(png, dpi=60, bbox_inches="tight")
