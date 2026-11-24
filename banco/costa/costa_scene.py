"""Alcuni passaggi singoli: il vento del radar, di ECMWF e di Lens, in nodi (per vedere a occhio scie e accelerazioni).
  python3 costa_scene.py <dati> <out> <png> area settore quante"""
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_analisi as ca
import costa_figura as cf

KN = 1.943844


def main():
    dati, out, png, area, sett, quante = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lat, lon, scene, _ = ca.carica(dati, area)
    ny, nx = len(lat), len(lon)
    R = np.array([np.hypot(s["su"], s["sv"]) for s in scene]); M = np.array([np.hypot(s["mu"], s["mv"]) for s in scene])
    raw = np.fromfile(os.path.join(out, f"lens_{area}.bin"), dtype="<f4", offset=16).reshape(len(scene), 2, ny, nx)
    L = np.hypot(raw[:, 0], raw[:, 1])
    mai = ~np.isfinite(R).any(axis=0)
    mu = np.array([np.nanmean(s["mu"]) for s in scene]); mv = np.array([np.nanmean(s["mv"]) for s in scene])
    da = (np.degrees(np.arctan2(-mu, -mv)) + 360) % 360
    q = cf.SETT.index(sett)
    sel = np.where((((da + 22.5) // 45).astype(int) % 8 == q) & (np.isfinite(R).sum(axis=(1, 2)) > 0.8 * (~mai).sum()))[0]
    sel = sel[np.argsort(-np.nanmean(M[sel], axis=(1, 2)))][:quante]
    fig, ax = plt.subplots(len(sel), 3, figsize=(10, 3.1 * len(sel) * np.ptp(lat) / np.ptp(lon) / np.cos(np.radians(lat.mean())) + 0.8), squeeze=False)
    fig.patch.set_facecolor("#fcfcfb")
    ext = [lon[0], lon[-1], lat[0], lat[-1]]; asp = 1 / np.cos(np.radians(lat.mean()))
    for r, k in enumerate(sel):
        vmax = np.nanpercentile(R[k] * KN, 99)
        for c, (F, nome) in enumerate([(R[k], "radar"), (M[k], "ECMWF"), (np.where(mai, np.nan, L[k]), "Lens")]):
            a = ax[r][c]; a.set_facecolor("#d9d7d0")
            im = a.imshow(np.ma.masked_invalid(F * KN), origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=vmax, aspect=asp, interpolation="nearest")
            a.contour(lon, lat, mai.astype(float), [0.5], colors="white", linewidths=0.5)
            a.set_xticks([]); a.set_yticks([])
            t = np.datetime64(int(scene[k]["t"]), "s")
            a.set_title(f"{nome} · {str(t)[:16].replace('T', ' ')} UTC", fontsize=8)
        fig.colorbar(im, ax=ax[r].tolist(), shrink=0.8, pad=0.01, label="nodi")
    fig.savefig(os.path.join(png, f"scene_{area}_{sett}.png"), dpi=110, bbox_inches="tight", facecolor=fig.get_facecolor())


if __name__ == "__main__":
    main()
