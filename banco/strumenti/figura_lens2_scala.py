"""Come si vede il primo Lens a 100 m con Lens 2 sopra: niente (la terza colonna), relativo per cella (Shift 0.58),
e a grande scala (il radar fissa solo il livello oltre s km). Vento uniforme, fattore della libreria per la stessa
situazione. Uso:
  python3 figura_lens2_scala.py <campi fini .bin> <indice vento> <atlante fine .lna> <libreria .lra>
                                <da (gradi)> <nodi> <v925 m/s> <scala km> <png> [S N O E della vista]
"""
import os
import sys

import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

qui = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [qui, os.path.join(qui, "..", "cfd")]
from controlla_atlante import leggi as leggi_lna   # noqa: E402
from leggi_campi import leggi                       # noqa: E402
from libreria_radar import Libreria                 # noqa: E402

p_bin, iv, p_lna, p_lib, da, kn, v925, scala, png = sys.argv[1:10]
vista = [float(x) for x in sys.argv[10:14]] if len(sys.argv) >= 14 else None
iv, da, kn, v925, scala = int(iv), float(da), float(kn), float(v925), float(scala)
la, lo, F = leggi(p_bin, iv + 1); U, V = F[iv][1]
A = leggi_lna(p_lna); m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
S = np.where(mare, np.hypot(U, V) * 1.944, np.nan)                     # nodi
L = Libreria(p_lib)
a = np.radians(da); s0 = kn * 0.5144
g = L.fattore(s0 * np.sin(a), s0 * np.cos(a), v925 * np.sin(a), v925 * np.cos(a), 0.0, np.full(L.nc, da))
G = np.full((L.ny, L.nx), np.nan); G.flat[np.array(L.celle)] = g
rlat = L.lat0 + np.arange(L.ny) * L.dlat; rlon = L.lon0 + np.arange(L.nx) * L.dlon
# cella del radar di ogni punto fine
jr = np.clip(np.rint((la - L.lat0) / L.dlat).astype(int), 0, L.ny - 1); ir = np.clip(np.rint((lo - L.lon0) / L.dlon).astype(int), 0, L.nx - 1)
JR, IR = np.meshgrid(jr, ir, indexing="ij")
# 1) relativo per cella: g / (media di Lens nella cella / modello), tra 0,5 e 2; dove la libreria non ha la cella, 1
R = S / kn
somma = np.zeros((L.ny, L.nx)); n = np.zeros((L.ny, L.nx))
ok = np.isfinite(R)
np.add.at(somma, (JR[ok], IR[ok]), R[ok]); np.add.at(n, (JR[ok], IR[ok]), 1)
rcella = np.where(n >= 5, somma / np.maximum(n, 1), np.nan)
g2 = np.clip(G / rcella, 0.5, 2.0); g2 = np.where(np.isfinite(g2), g2, 1.0)
S_rel = S * g2[JR, IR]
# 2) grande scala: exp(liscio(log g)) / liscio(Lens / modello), sigma = scala / 2
dkm = L.dlat * 111.19


def liscia(Z, sig):
    mk = np.isfinite(Z); z0 = np.where(mk, Z, 0.0)
    return gaussian_filter(z0, sig) / np.maximum(gaussian_filter(mk.astype(float), sig), 1e-6)


Gs = np.exp(liscia(np.log(G), scala / 2 / dkm))
fine_km = abs(la[1] - la[0]) * 111.19
Rs = liscia(R, scala / 2 / fine_km)
S_sca = S * Gs[JR, IR] / Rs
fig, ax = plt.subplots(1, 3, figsize=(24, 8))
titoli = ["primo Lens a 100 m, da solo (la terza colonna)", "con Lens 2 relativo cella per cella (Shift 0.58)",
          f"con Lens 2 a grande scala ({scala:g} km): il radar fissa solo il livello"]
for k, (a_, Z, t) in enumerate(zip(ax, (S, S_rel, S_sca), titoli)):
    im = a_.pcolormesh(lo, la, Z, cmap="turbo", vmin=0, vmax=24, shading="auto", rasterized=True)
    a_.contour(lo, la, mare.astype(float), [0.5], colors="k", linewidths=0.5)
    a_.set_title(f"da {da:g}°, {kn:g} nodi — {t}", fontsize=10); a_.set_aspect(1 / np.cos(np.radians(la.mean())))
    if vista:
        a_.set_ylim(vista[0], vista[1]); a_.set_xlim(vista[2], vista[3])
plt.colorbar(im, ax=list(ax), shrink=0.7, label="nodi", pad=0.01)
plt.savefig(png, dpi=60, bbox_inches="tight")
for nome, Z in (("primo", S), ("relativo", S_rel), ("scala", S_sca)):
    print(f"{nome:9s} mediana {np.nanmedian(Z):.1f}  quota sotto 10 nodi {100 * np.nanmean(Z < 10):.0f}%  sotto 6 {100 * np.nanmean(Z < 6):.0f}%")
