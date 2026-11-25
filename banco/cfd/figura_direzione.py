"""Figura: la carta di Pietro sotto il Pellegrino com'è e con la direzione del rilievo pieno (LENS_DIR_RILIEVO=1),
forza invariata. La rotazione è quella dell'atlante di Lens per il vento dell'ora (classe Fr 0,7)."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "banco/strumenti"); sys.path.insert(0, "banco/cfd")
from controlla_atlante import leggi
from taratura_radar import mare_su
atl, binp, ora, uscita = sys.argv[1:5]
ora = int(ora)
A = leggi(atl); s, n, w, e = A["box"]; ny, nx = A["ny"], A["nx"]
lat = n - (np.arange(ny) + 0.5) * (n - s) / ny; lon = w + (np.arange(nx) + 0.5) * (e - w) / nx
mare = mare_su(atl, lat, lon)
b = np.load(binp); V = b[0, ora].reshape(ny, nx) / 2.0; D = b[2, ora].reshape(ny, nx) * 360 / 128
d0 = float(np.degrees(np.arctan2(np.mean(np.sin(np.radians(D[mare]))), np.mean(np.cos(np.radians(D[mare]))))) % 360)
k = int(np.argmin(np.abs(np.array(A["fr"]) - 0.7))); ue, ve, un, vn = A["campi"][k]
u0, v0 = -np.sin(np.radians(d0)), -np.cos(np.radians(d0))
U = u0 * ue + v0 * un; Vv = u0 * ve + v0 * vn
rot = (np.degrees(np.arctan2(-U, -Vv)) - d0 + 540) % 360 - 180
D2 = (D + rot) % 360
fig, ax = plt.subplots(1, 2, figsize=(15, 7.5))
LA, LO = np.meshgrid(lat, lon, indexing="ij")
sel = np.zeros_like(mare); sel[::3, ::3] = True; sel &= mare
for a, (tit, DD) in zip(ax, [(f"Com'è oggi (vento al largo da {d0:.0f}°)", D), ("Con la direzione del rilievo pieno, stessa forza", D2)]):
    im = a.pcolormesh(lon, lat, np.where(mare, V, np.nan), cmap="viridis", vmin=4, vmax=16, shading="auto")
    a.contour(lon, lat, mare.astype(float), [0.5], colors="k", linewidths=0.7)
    a.quiver(LO[sel], LA[sel], -np.sin(np.radians(DD[sel])), -np.cos(np.radians(DD[sel])), scale=40, width=0.0025)
    a.set_xlim(13.33, 13.43); a.set_ylim(38.13, 38.23); a.set_aspect(1 / np.cos(np.radians(38.18))); a.set_title(tit, fontsize=10)
    for nome, (la_, lo_) in {"P. del Rotolo": (38.192, 13.362), "Vergine Maria": (38.170, 13.368), "Arenella": (38.150, 13.373)}.items():
        a.plot(lo_, la_, "r.", ms=5); a.text(lo_ + 0.002, la_, nome, fontsize=7, color="r")
plt.colorbar(im, ax=ax, shrink=0.7, label="nodi")
plt.savefig(uscita, dpi=90, bbox_inches="tight")
for nome, (la_, lo_) in {"Punta del Rotolo, 400 m": (38.193, 13.367), "Vergine Maria, 400 m": (38.170, 13.373), "Arenella, 400 m": (38.150, 13.380)}.items():
    j = int(round((n - la_) / (n - s) * ny - 0.5)); i = int(round((lo_ - w) / (e - w) * nx - 0.5))
    print(f"{nome:24s} {V[j, i]:.1f} kn, da {D[j, i]:.0f}° -> {D2[j, i]:.0f}°")
