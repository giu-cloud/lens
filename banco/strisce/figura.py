"""Figura di un passaggio: l'immagine radar (senza il fondo lento) con le direzioni lette dalle strisce e quelle del modello."""
import sys, os, glob, datetime as dt, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zoneinfo import ZoneInfo
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analizza import immagine, terra, D
from strisce import prepara, liscia

area, chiave, usc = sys.argv[1], sys.argv[2], sys.argv[3]
QMIN = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
M = terra(area)
f = [x for x in sorted(glob.glob(f"{D}/strisce/{area}_*.npz")) if chiave in x][0]
x, g = immagine(area, f, M)
b, fz = prepara(g["amp"], cella_m=200.0); b = liscia(b, 0.8)
ny, nx = b.shape
fig, ax = plt.subplots(figsize=(15, 15 * ny / nx + 1.2))
v = np.nanstd(b)
ax.imshow(b, cmap="gray", vmin=1 - 3 * v, vmax=1 + 3 * v, extent=(0, nx * fz, ny * fz, 0), interpolation="nearest")
ax.set_facecolor("#E9E4D8")
r, c, d = g["d"]["r"], g["d"]["c"], g["d"]
passo = 4
sel = np.zeros(r.shape, bool); sel[passo // 2::passo, passo // 2::passo] = True
ok = sel & (d["n"] > 0.7) & (d["q"] >= QMIN)
a = np.radians(d["ang"]); L = 45
ax.quiver(c[ok], r[ok], np.cos(a[ok]) * L, np.sin(a[ok]) * L, angles="xy", scale_units="xy", scale=1, headwidth=0, headlength=0, headaxislength=0, pivot="mid", color="#FFB000", width=0.0022, label="strisce del radar (qualità sufficiente)")
# il modello nel piano dell'immagine: si risolve [ec er] (cos a, sin a) = (sen b, cos b), b = direzione di bussola
def nel_piano(bussola):
    b = np.radians(bussola); tx, ty = np.sin(b), np.cos(b)
    ec, er = g["ec"], g["er"]; det = ec[..., 0] * er[..., 1] - er[..., 0] * ec[..., 1]
    return np.arctan2((ec[..., 0] * ty - ec[..., 1] * tx) / det, (tx * er[..., 1] - ty * er[..., 0]) / det)
da = nel_piano(g["dm"])
ok2 = sel & (d["n"] > 0.7)
ax.quiver(c[ok2], r[ok2], np.cos(da[ok2]) * L, np.sin(da[ok2]) * L, angles="xy", scale_units="xy", scale=1, headwidth=0, headlength=0, headaxislength=0, pivot="mid", color="#2B8CFF", width=0.0012, alpha=0.8, label="ECMWF")
ax.set_xticks([]); ax.set_yticks([])
t = dt.datetime.fromtimestamp(g["t"], dt.timezone.utc)
ax.set_title(f"Sentinel-1, {t.astimezone(ZoneInfo('Europe/Rome')):%d/%m/%Y %H:%M} ora italiana. Immagine radar a celle di 200 m, tolto il fondo lento.\nTrattini gialli: direzione letta dalle strisce. Trattini azzurri: direzione di ECMWF. (Immagine nel verso del satellite, non a nord.)", fontsize=11)
ax.legend(loc="lower right", fontsize=9)
fig.savefig(usc, dpi=80, bbox_inches="tight"); print(usc)
