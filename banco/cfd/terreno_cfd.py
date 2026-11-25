"""Il terreno per la CFD di un sito, indipendente dagli atlanti di Lens (che coprono solo il riquadro dell'area).

  python3 terreno_cfd.py <uscita .npz> --centro LAT,LON [--raggio 26] [--dx 100]

Quadrato di lato 2 x raggio km attorno al centro: quota da Copernicus GLO-30 (media), copertura del suolo da ESA
WorldCover (moda) e la stessa scabrezza z0 di Lens (lens/calcola.py). Righe da nord. Il raggio deve coprire gli angoli
del dominio ruotato della CFD (per un dominio di +-18 km servono almeno 18 x 1,42 = 25,5 km).
"""
import argparse
import math
import os
import sys

import numpy as np

os.environ.setdefault("LENS_PAGINA", "0")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from lens import calcola as C   # noqa: E402

Z0 = {10: 0.6, 20: 0.15, 30: 0.05, 40: 0.08, 50: 0.8, 60: 0.02, 80: 0.0005, 90: 0.05}


def main():
    from rasterio.enums import Resampling
    ap = argparse.ArgumentParser()
    ap.add_argument("uscita"); ap.add_argument("--centro", required=True)
    ap.add_argument("--raggio", type=float, default=26.0); ap.add_argument("--dx", type=float, default=100.0)
    a = ap.parse_args()
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    kx = 111.19 * math.cos(math.radians(lat0))
    S, N = lat0 - a.raggio / 111.19, lat0 + a.raggio / 111.19
    O, E = lon0 - a.raggio / kx, lon0 + a.raggio / kx
    ny = nx = int(round(2 * a.raggio * 1000 / a.dx))
    box = (S, N, O, E)
    h = np.maximum(C.mosaico(C.url_dem, 1, Resampling.average, np.float32, box, (ny, nx)), 0.0)
    try:
        c = C.mosaico(C.url_cov, 3, Resampling.mode, np.uint8, box, (ny, nx))
    except Exception as ex:
        print("WorldCover non letto:", ex)
        c = np.zeros((ny, nx), np.uint8)
    mare = (h <= 0.5) & ((c == 80) | (c == 0))
    z0 = np.full(h.shape, 0.10, np.float32)
    for classe, v in Z0.items():
        z0[c == classe] = v
    z0[mare] = 0.0005
    lat = N - (np.arange(ny) + 0.5) * (N - S) / ny
    lon = O + (np.arange(nx) + 0.5) * (E - O) / nx
    np.savez_compressed(a.uscita, h=h.astype("f4"), z0=z0, mare=mare, lat=lat, lon=lon, centro=np.array([lat0, lon0]))
    print(f"{a.uscita}: {ny} x {nx} celle da {a.dx:.0f} m, {S:.3f}-{N:.3f} N, {O:.3f}-{E:.3f} E, "
          f"quota massima {h.max():.0f} m, mare {mare.mean() * 100:.0f}%")


if __name__ == "__main__":
    main()
