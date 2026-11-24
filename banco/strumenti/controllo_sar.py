"""Uso: python3 controllo_sar.py <cartella con i costa_palermo_*.npz della release costa-sar di shift-grib>

Controllo rapido: rapporto mediano radar / ECMWF (quello del prodotto) per distanza dalla terra, Palermo."""
import sys, glob, os
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "costa"))
import costa_analisi as ca
lat, lon, scene, info = ca.carica(sys.argv[1], "palermo")
print(len(scene), "scene")
S = np.array([np.hypot(s["su"], s["sv"]) for s in scene]); M = np.array([np.hypot(s["mu"], s["mv"]) for s in scene])
mare = np.isfinite(S).any(axis=0)               # celle mai misurate = terra (approssimazione)
dy = 111.0 * abs(lat[1] - lat[0]); dx = 111.32 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
dist = ndimage.distance_transform_edt(mare, sampling=(dy, dx))
print("celle mare", mare.sum(), "di", mare.size)
ok = np.isfinite(S) & np.isfinite(M) & (M * 1.943844 > 6)
for a, b in [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 99)]:
    sel = ok & (dist[None] > a) & (dist[None] <= b)
    if sel.sum() > 50:
        print(f"{a:4}-{b:<4} km  celle {sel.sum():6d}  radar/ECMWF mediano {np.median(S[sel] / M[sel]):.2f}")
ore = np.array([ (s["t"] // 3600) % 24 for s in scene]); print("ore UTC dei passaggi:", sorted(set(ore.tolist())))
