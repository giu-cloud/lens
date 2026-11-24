"""Come diagnosi_direzioni.py, ma con il vento ECMWF che accompagna ogni passaggio nel prodotto del radar (circa 0,1°,
lettura del prodotto): tutte le scene, non solo quelle con l'archivio di shift-grib. Serve a vedere i fenomeni con più
casi; non è la lettura di Shift.

  python3 diagnosi_direzioni_cmems.py <cartella dei costa_*.npz> area [area ...]
"""
import os, sys
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "costa"))
import costa_analisi as ca
KN = 1.943844
SETT = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO"]
for area in sys.argv[2:]:
    lat, lon, scene, _ = ca.carica(sys.argv[1], area)
    R = np.array([np.hypot(s["su"], s["sv"]) for s in scene]); M = np.array([np.hypot(s["mu"], s["mv"]) for s in scene])
    mu = np.array([s["mu"] for s in scene]); mv = np.array([s["mv"] for s in scene])
    mare = np.isfinite(R).any(axis=0)
    kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
    D = np.broadcast_to(ndimage.distance_transform_edt(mare, sampling=(kmy, kmx)), R.shape)
    larg = (D > 8) & np.isfinite(M)
    U = np.array([np.nanmedian(M[i][larg[i]]) if larg[i].any() else np.nan for i in range(len(R))])
    dd = np.array([(np.degrees(np.arctan2(-np.nanmedian(mu[i][larg[i]]), -np.nanmedian(mv[i][larg[i]]))) + 360) % 360 if larg[i].any() else np.nan for i in range(len(R))])
    k = np.where(np.isfinite(dd), np.round(np.nan_to_num(dd) / 22.5).astype(int) % 16, -1)
    print(f"\n== {area}: {len(scene)} scene; per direzione del modello al largo (oltre 8 km), scene con oltre 8 nodi al largo")
    print(f"   {'da':<4} {'scene':>5}  {'0-3 km':>14}  {'3-8 km':>14}  {'oltre 8':>14}   (radar/modello, errore nodi)")
    for s in range(16):
        sc = (k == s) & (U * KN > 8)
        if sc.sum() < 4:
            continue
        cols = []
        for lo, hi in [(0, 3), (3, 8), (8, 99)]:
            m = sc[:, None, None] & np.isfinite(R) & np.isfinite(M) & (M > 1.5) & (D >= lo) & (D < hi)
            cols.append(f"{np.median(R[m] / M[m]):5.2f} {np.mean(M[m] - R[m]) * KN:+5.1f}" if m.sum() > 50 else "      -     ")
        print(f"   {SETT[s]:<4} {sc.sum():5d}  " + "  ".join(f"{c:>14}" for c in cols))
