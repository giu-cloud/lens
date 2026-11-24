"""Figura riassuntiva dello studio del 24-25/11: errore sistematico entro 3 km dalla costa, nodi (vento vero)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURF, INK, INK2, BLU, GRIGIO = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6", "#a8a7a2"
pannelli = [
    ("Aree mai viste (un'area lasciata fuori alla volta)", [
        ("Modello letto dal mare", 1.738, True), ("Lens 0.54", 1.693, False),
        ("Vecchio Lens: risposta del rilievo", 2.040, False), ("Vecchio Lens: ridosso lungo", 2.013, False),
        ("Moduli fisici combinati", 1.728, False), ("Alberi: rilievo e scena", 1.664, False),
        ("Alberi: rilievo, scena, fisica", 1.642, False)]),
    ("Aree mai viste, settimane della cassaforte", [
        ("Modello letto dal mare", 1.626, True), ("Lens 0.54", 1.581, False),
        ("Alberi: rilievo, scena, stabilità, fisica", 1.460, False)]),
    ("Aree note, settimane della cassaforte", [
        ("Modello letto dal mare", 1.626, True), ("Lens 0.54", 1.581, False),
        ("Carta imparata dal radar (cella × ottante)", 0.937, False)]),
]
fig, axs = plt.subplots(3, 1, figsize=(9, 7.2), gridspec_kw={"height_ratios": [7, 3, 3]}, facecolor=SURF)
for ax, (tit, righe) in zip(axs, pannelli):
    ax.set_facecolor(SURF)
    y = list(range(len(righe)))[::-1]
    for yy, (nome, v, rif) in zip(y, righe):
        ax.barh(yy, v, height=0.62, color=GRIGIO if rif else BLU)
        ax.text(v + 0.025, yy, f"{v:.2f}".replace(".", ","), va="center", fontsize=9, color=INK, zorder=3, bbox=dict(facecolor=SURF, edgecolor="none", pad=1))
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in righe], fontsize=9, color=INK)
    ax.set_xlim(0, 2.3); ax.set_title(tit, fontsize=10, loc="left", color=INK)
    ax.axvline(righe[0][1], color=INK2, lw=0.8, ls=":", zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRIGIO); ax.tick_params(axis="x", colors=INK2, labelsize=8); ax.tick_params(axis="y", length=0)
axs[-1].set_xlabel("errore sistematico entro 3 km dalla costa (nodi, vento vero a 10 m; meno è meglio)", fontsize=9, color=INK2)
fig.suptitle("Quanto ci si avvicina al radar: 10 aree, 3.750 passaggi dei Sentinel-1, ECMWF letto dal mare", fontsize=11, color=INK, x=0.02, ha="left")
plt.tight_layout()
plt.savefig("/home/claude/lens/banco/figure/metodi_a_confronto.png", dpi=110, facecolor=SURF)
