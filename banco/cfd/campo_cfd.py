"""Il vento a 10 m di un caso CFD (caso_cfd.py) sulla griglia interna dell'atlante di Lens, in rapporto al vento d'ingresso.

  python3 campo_cfd.py <atlante .lna o terreno .npz> <cartella del caso> <uscita .npz> --centro LAT,LON --da GRADI
         [--raggio-uscita 18 --dx 150]   (con il terreno .npz di un sito: griglia regolare attorno al centro)

Scrive u, v (est, nord) per 1 m/s di vento a 10 m sul mare all'ingresso, sulla griglia dell'atlante o del sito (righe da nord),
NaN fuori dal dominio della CFD o entro 2 km dai suoi bordi (dove il rilievo è sfumato).
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
from caso_cfd import Geometria
from confronto_cfd import cfd
from controlla_atlante import leggi as leggi_lna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlante"); ap.add_argument("caso"); ap.add_argument("uscita")
    ap.add_argument("--centro", required=True); ap.add_argument("--da", type=float, required=True)
    ap.add_argument("--raggio-uscita", type=float, default=18.0); ap.add_argument("--dx", type=float, default=150.0)
    a = ap.parse_args()
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    G = Geometria(lat0, lon0, a.da)
    bm = open(f"{a.caso}/system/blockMeshDict").read()
    nx, ny, nz = (int(v) for v in re.search(r"hex \(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)", bm).groups())
    vin = float(re.search(r"Uref ([\d.]+)", open(f"{a.caso}/0/U").read()).group(1))
    xr, yr, U10, it = cfd(a.caso, nx, ny, nz)
    if a.atlante.endswith(".npz"):
        # sito senza atlante: griglia regolare attorno al centro, lato 2 x raggio, con la maglia della CFD
        dxm = a.dx
        R = a.raggio_uscita; kx = 111.19 * np.cos(np.radians(lat0)); nn = int(round(2 * R * 1000 / dxm))
        lat = lat0 + R / 111.19 - (np.arange(nn) + 0.5) * (2 * R / 111.19) / nn
        lon = lon0 - R / kx + (np.arange(nn) + 0.5) * (2 * R / kx) / nn
    else:
        L = leggi_lna(a.atlante)
        s, n, w, e = L["box"]
        lat = n - (np.arange(L["ny"]) + 0.5) * (n - s) / L["ny"]; lon = w + (np.arange(L["nx"]) + 0.5) * (e - w) / L["nx"]
    LO, LA = np.meshgrid(lon, lat)
    xp, yp = G.ruotato(LA, LO)
    dx = xr[0, 1] - xr[0, 0]; dy = yr[1, 0] - yr[0, 0]
    fi = (xp - xr[0, 0]) / dx; fj = (yp - yr[0, 0]) / dy
    dentro = (xp > xr.min() + 2) & (xp < xr.max() - 2) & (yp > yr.min() + 2) & (yp < yr.max() - 2)
    i0 = np.clip(np.floor(fi).astype(int), 0, nx - 2); j0 = np.clip(np.floor(fj).astype(int), 0, ny - 2)
    tx = np.clip(fi - i0, 0, 1)[..., None]; ty = np.clip(fj - j0, 0, 1)[..., None]
    Ub = ((1 - ty) * (1 - tx) * U10[j0, i0] + (1 - ty) * tx * U10[j0, i0 + 1] + ty * (1 - tx) * U10[j0 + 1, i0] + ty * tx * U10[j0 + 1, i0 + 1]) / vin
    u = Ub[..., 0] * G.f[0] + Ub[..., 1] * G.l[0]; v = Ub[..., 0] * G.f[1] + Ub[..., 1] * G.l[1]
    u[~dentro] = np.nan; v[~dentro] = np.nan
    np.savez_compressed(a.uscita, u=u.astype("f4"), v=v.astype("f4"), da=a.da, lat=lat, lon=lon, iterazione=int(it))
    print(f"{a.uscita}: direzione {a.da:.1f}°, {np.isfinite(u).mean() * 100:.0f}% della griglia, iterazione {it}")


if __name__ == "__main__":
    main()
