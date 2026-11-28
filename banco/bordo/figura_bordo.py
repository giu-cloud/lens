"""Figura del bordo: Modelli | Lens di prima (0.59.1) | Lens con l'alone (0.59.2), più il profilo lungo una linea che
attraversa il bordo. Dati da FiguraBordo.kt, colori del vento dell'app."""
import sys, datetime as dt, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from zoneinfo import ZoneInfo

dati, terra, usc, t, titolo, lon_p = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5], float(sys.argv[6])
STOP = [(0, '#7FB2E5'), (5, '#4CAF50'), (10, '#C0CA33'), (15, '#FFB300'), (20, '#FF7043'), (25, '#E53935'), (32, '#8E24AA')]
VENTO = LinearSegmentedColormap.from_list('app', [(k / 32, c) for k, c in STOP])
KN = 1.943844
b = open(dati, "rb").read(); assert b[:4] == b"FBD1"
NY, NX = np.frombuffer(b, "<i4", 2, 4); la0, lo0, ps = np.frombuffer(b, "<f8", 3, 12); box = np.frombuffer(b, "<f8", 4, 36)
F = np.frombuffer(b, "<f4", 6 * NY * NX, 68).reshape(6, NY, NX)
lat = la0 + np.arange(NY) * ps; lon = lo0 + np.arange(NX) * ps
rb = open(terra, "rb").read(); ny, nx = np.frombuffer(rb, ">i4", 2); M = np.frombuffer(rb, np.uint8, offset=8).reshape(ny, nx)
mla = 37.70 + (np.arange(ny) + 0.5) * 0.002; mlo = 15.20 + (np.arange(nx) + 0.5) * 0.002
J = np.clip(np.rint((lat - mla[0]) / 0.002).astype(int), 0, ny - 1); I = np.clip(np.rint((lon - mlo[0]) / 0.002).astype(int), 0, nx - 1)
terra_g = M[J][:, I].astype(bool)
LA, LO = np.meshgrid(lat, lon, indexing="ij")
d = dt.datetime.fromtimestamp(t, dt.timezone.utc)
fig = plt.figure(figsize=(26, 8.6)); gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15], wspace=0.08)
nomi = ["Modelli (ECMWF)", "Lens di prima (0.59.1): Lens dentro, il modello fuori", "Lens con l'alone (0.59.2)"]
for k in range(3):
    ax = fig.add_subplot(gs[k]); u, v = F[2 * k], F[2 * k + 1]; s = np.hypot(u, v) * KN
    ax.pcolormesh(lon, lat, np.where(terra_g, np.nan, s), cmap=VENTO, vmin=0, vmax=32, shading="auto", rasterized=True)
    ax.contourf(lon, lat, terra_g.astype(float), [0.5, 1.5], colors=["#F3F0E8"])
    ax.contour(lon, lat, terra_g.astype(float), [0.5], colors="#8EA2B2", linewidths=0.7)
    q = np.zeros_like(terra_g); q[6::12, 6::12] = True; q &= ~terra_g
    sp = np.hypot(u, v) + 1e-9
    ax.quiver(LO[q], LA[q], (u / sp)[q], (v / sp)[q], scale=40, width=0.004, color="#1F2D3A", alpha=0.65)
    ax.plot([box[2], box[3], box[3], box[2], box[2]], [box[0], box[0], box[1], box[1], box[0]], "--", color="#1F2D3A", lw=1.2)
    ax.axvline(lon_p, color="white", lw=1.0, alpha=0.9)
    ax.set_aspect(1 / np.cos(np.radians(lat.mean()))); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(lon[0], lon[-1]); ax.set_ylim(lat[0], lat[-1])
    ax.set_title(nomi[k], fontsize=12)
ax = fig.add_subplot(gs[3])
i = int(np.argmin(np.abs(lon - lon_p)))
km = (lat - box[0]) * 111.19
for k, (nm, c, st) in enumerate(zip(["Modelli", "Lens di prima", "Lens con l'alone"], ["#7F8C8D", "#E67E22", "#2471A3"], ["-", "--", "-"])):
    s = np.hypot(F[2 * k][:, i], F[2 * k + 1][:, i]) * KN
    s[terra_g[:, i]] = np.nan
    ax.plot(km, s, st, color=c, lw=2.2, label=nm)
ax.axvline(0, color="#1F2D3A", ls="--", lw=1); ax.axvspan(0, 1.5 * 1.852, color="#1F2D3A", alpha=0.06)
ax.text(0.3, ax.get_ylim()[0] + 0.3, "bordo del riquadro;\nin grigio il raccordo\ndi 1,5 miglia", fontsize=10, va="bottom")
ax.set_xlim(-16, 16); ax.set_xlabel("km dal bordo sud del riquadro (negativo = fuori)"); ax.set_ylabel("nodi")
ax.set_title(f"Lungo la linea bianca ({lon_p:.2f} E)", fontsize=12); ax.legend(loc="upper left"); ax.grid(alpha=0.3)
fig.suptitle(f"{titolo}, {d.astimezone(ZoneInfo('Europe/Rome')):%d/%m/%Y %H:%M} ora italiana ({d:%H:%M} UTC). Previsione ECMWF delle 12 UTC; tratteggio = riquadro di Lens.", fontsize=14, y=0.99)
fig.savefig(usc, dpi=62, bbox_inches="tight"); print(usc)
