"""Le mappe dell'atlante: per direzione del vento, quanto vento in più o in meno del modello vede il radar sotto costa,
e che cosa fa Lens. Usa i campi Lens già calcolati dalla diagnosi (out/lens_<area>.bin).

  python3 costa_figura.py <dati> <out> <cartella png> area [area ...]
"""
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_analisi as ca

KN = 1.943844
SETT = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
TITOLI = {"palermo": "Golfo di Palermo", "egadi": "Egadi", "eolie": "Eolie", "napoli": "Golfo di Napoli", "toscano": "Arcipelago Toscano",
          "maddalena-bonifacio": "La Maddalena e Bonifacio", "messina": "Stretto di Messina", "malta": "Malta", "cagliari": "Golfo di Cagliari", "genova": "Genova"}


def mappa_colori():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("div", ["#123f73", "#2a78d6", "#9cc0ea", "#f0efec", "#f0a3a0", "#e34948", "#8a1f1f"])


def main():
    dati, out, png = sys.argv[1:4]
    os.makedirs(png, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#b9b8b3", "axes.linewidth": 0.6})
    for area in sys.argv[4:]:
        lat, lon, scene, _ = ca.carica(dati, area)
        ny, nx = len(lat), len(lon)
        SU = np.array([s["su"] for s in scene]); SV = np.array([s["sv"] for s in scene])
        MU = np.array([s["mu"] for s in scene]); MV = np.array([s["mv"] for s in scene])
        R = np.hypot(SU, SV); M = np.hypot(MU, MV)
        lb = os.path.join(out, f"lens_{area}.bin")
        L = None
        if os.path.exists(lb):
            raw = np.fromfile(lb, dtype="<f4", offset=16)
            if raw.size == len(scene) * 2 * ny * nx:
                raw = raw.reshape(len(scene), 2, ny, nx); L = np.hypot(raw[:, 0], raw[:, 1]).astype("f8")
        mai = ~np.isfinite(R).any(axis=0)
        kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
        dist = ndimage.distance_transform_edt(~mai, sampling=(kmy, kmx))
        # direzione e forza del modello su tutta l'area, per scena
        mu = np.nanmean(MU, axis=(1, 2)); mv = np.nanmean(MV, axis=(1, 2)); ms = np.nanmean(M, axis=(1, 2))
        da = (np.degrees(np.arctan2(-mu, -mv)) + 360) % 360
        sett = ((da + 22.5) // 45).astype(int) % 8
        righe = [("radar", "Radar ÷ ECMWF\n(misurato)")] + ([("lens", "Lens ÷ ECMWF\n(calcolato)")] if L is not None else [])
        usati = [q for q in range(8) if ((sett == q) & (ms > 2.5)).sum() >= 6]
        print(area, "passaggi per settore (oltre 5 nodi):", {SETT[q]: int(((sett == q) & (ms > 2.5)).sum()) for q in range(8)})
        if not usati:
            continue
        nc = min(4, len(usati)); gruppi = (len(usati) + nc - 1) // nc; nr = len(righe)
        alto = np.ptp(lat) / np.ptp(lon) / np.cos(np.radians(lat.mean()))
        fig, ax = plt.subplots(nr * gruppi, nc, figsize=(2.6 * nc + 1.6, 2.6 * alto * nr * gruppi + 0.5 * gruppi + 0.8), squeeze=False)
        fig.patch.set_facecolor("#fcfcfb")
        ext = [lon[0], lon[-1], lat[0], lat[-1]]
        asp = 1 / np.cos(np.radians(lat.mean()))
        norm = TwoSlopeNorm(1.0, 0.5, 1.5)
        for a in ax.flat:
            a.set_visible(False)
        for c0, q in enumerate(usati):
            g, c = divmod(c0, nc)
            idx = np.where((sett == q) & (ms > 2.5))[0]
            for r, (k, nome) in enumerate(righe):
                F = R if k == "radar" else L
                ok = np.isfinite(F[idx]) & np.isfinite(M[idx]) & (M[idx] > 2.0)
                num = np.where(ok, F[idx], 0).sum(axis=0); den = np.where(ok, M[idx], 0).sum(axis=0); cnt = ok.sum(axis=0)
                rap = np.where((cnt >= 4) & ~mai, num / np.maximum(den, 1e-6), np.nan)
                a = ax[g * nr + r][c]
                a.set_visible(True)
                a.set_facecolor("#d9d7d0")          # terra e celle senza dati
                im = a.imshow(np.ma.masked_invalid(rap), origin="lower", extent=ext, cmap=mappa_colori(), norm=norm, aspect=asp, interpolation="nearest")
                a.contour(lon, lat, mai.astype(float), [0.5], colors="#3d3d3a", linewidths=0.5)
                a.set_xticks([]); a.set_yticks([])
                if r == 0:
                    a.set_title(f"vento da {SETT[q]} · {len(idx)} passaggi · {np.mean(ms[idx]) * KN:.0f} nodi", fontsize=8.5, color="#3d3d3a")
                if c == 0:
                    a.set_ylabel(nome, fontsize=8.5, color="#3d3d3a")
        fig.subplots_adjust(left=0.07, right=0.86, top=1 - 0.75 / fig.get_figheight(), bottom=0.02, wspace=0.05, hspace=0.12)
        cax = fig.add_axes([0.885, 0.2, 0.018, 0.6])
        cb = fig.colorbar(im, cax=cax, ticks=[0.5, 0.75, 1.0, 1.25, 1.5])
        cb.ax.set_yticklabels(["½ del modello", "−25 %", "come il modello", "+25 %", "+50 %"], fontsize=7.5)
        cb.outline.set_linewidth(0.4)
        fig.suptitle(f"{TITOLI.get(area, area)}: vento sotto costa rispetto a ECMWF, per direzione (Sentinel-1, {len(scene)} passaggi 2024-2026)",
                     fontsize=10, color="#1f1f1d", x=0.02, ha="left", y=1 - 0.2 / fig.get_figheight())
        fig.savefig(os.path.join(png, f"atlante_{area}.png"), dpi=130, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)
        print(area, "fatto", [SETT[q] for q in usati])


if __name__ == "__main__":
    main()
