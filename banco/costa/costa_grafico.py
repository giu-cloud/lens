"""Un grafico: errore medio del vento (modello meno radar) per distanza dalla terra, tutte le aree insieme."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_pool as cp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nuova = sys.argv[1]           # variante di Lens scelta (cartella in outT)
a = cp.pool("outT/lens", ["shift", "mare", "lens", "lens_mare"])
b = cp.pool(f"outT/{nuova}", ["lens_mare"])
serie = [("Shift fino alla 0.53", a["shift"], "#2a78d6", "-", "o"),
         ("Shift 0.54: solo punti di mare", a["mare"], "#1baf7a", "-", "s"),
         ("Lens fino alla 0.53", a["lens"], "#eb6834", "--", "^"),
         ("Lens 0.54 (ritarato col radar)", b["lens_mare"], "#eda100", "--", "D")]
x = np.arange(5); etichette = ["0-1,5", "1,5-3", "3-5", "5-10", "oltre 10"]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
fig.patch.set_facecolor("#fcfcfb"); ax.set_facecolor("#fcfcfb")
ax.axhline(0, color="#8a8984", lw=0.8)
for nome, v, col, ls, mk in serie:
    y = [e for e, q, n in v]
    ax.plot(x, y, color=col, lw=2, ls=ls, marker=mk, ms=7, mec="#fcfcfb", mew=1.2, label=nome)
    ax.annotate(f"{y[0]:+.1f}".replace(".", ",").replace("-", "−"), (x[0], y[0]), xytext=(-30, 0), textcoords="offset points", va="center", fontsize=8.5, color="#3d3d3a")
ax.set_xticks(x, etichette); ax.set_xlabel("distanza dalla terra (km)", color="#3d3d3a")
ax.set_ylabel("errore medio del vento (nodi)\nnegativo = meno del vero", color="#3d3d3a")
ax.set_xlim(-0.6, 4.3)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("#b9b8b3")
ax.tick_params(colors="#3d3d3a"); ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
from matplotlib.ticker import FuncFormatter
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.1f}".replace(".", ",").replace("-", "−")))
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
n = sum(q for _, _, q in a["shift"])
ax.set_title(f"ECMWF letto da Shift contro il radar dei Sentinel-1\n10 aree, 1.054 passaggi, {n / 1e6:.2f} milioni di celle", fontsize=10, color="#1f1f1d", loc="left")
fig.savefig("png/errore_per_distanza.png", dpi=130, bbox_inches="tight", facecolor=fig.get_facecolor())
