"""La forza del vento di Lens e della CFD contro il radar, sulle celle che il radar vede (almeno 1,5-2 km dalla riva).

  python3 radar_contro_cfd.py <composito radar .npz> <atlante .lna> <cartella del caso> --centro LAT,LON --da GRADI [--fr 1.0]

Il composito è la media di radar ÷ vento al largo sui passaggi con il vento da quella direzione (radar_ponente_*.npz).
Si confronta il disegno, non il livello: ogni campo si divide per la sua mediana sulle celle comuni, poi correlazione
e scarto quadratico; anche solo nella zona sottovento al Pellegrino (a est di 13,36 E, a sud di 38,20 N).
"""
import argparse
import math
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
    ap.add_argument("radar"); ap.add_argument("atlante"); ap.add_argument("caso")
    ap.add_argument("--centro", required=True); ap.add_argument("--da", type=float, required=True); ap.add_argument("--fr", type=float, default=1.0)
    ap.add_argument("--vento", type=float, default=8.0)
    a = ap.parse_args()
    R = np.load(a.radar); lat, lon, rap = R["lat"], R["lon"], R["rapporto"]
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    ok = np.isfinite(rap)
    la_p, lo_p, r_p = LA[ok], LO[ok], rap[ok]
    # CFD
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    G = Geometria(lat0, lon0, a.da)
    bm = open(f"{a.caso}/system/blockMeshDict").read()
    nx, ny, nz = (int(v) for v in re.search(r"hex \(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)", bm).groups())
    xr, yr, U10, _ = cfd(a.caso, nx, ny, nz)
    sp = np.hypot(U10[..., 0], U10[..., 1]) / a.vento
    xp, yp = G.ruotato(la_p, lo_p)
    dentro = (xp > xr.min()) & (xp < xr.max()) & (yp > yr.min()) & (yp < yr.max())
    ii = np.clip(np.round((xp - xr[0, 0]) / (xr[0, 1] - xr[0, 0])).astype(int), 0, nx - 1)
    jj = np.clip(np.round((yp - yr[0, 0]) / (yr[1, 0] - yr[0, 0])).astype(int), 0, ny - 1)
    c_p = np.where(dentro, sp[jj, ii], np.nan)
    # Lens
    L = leggi_lna(a.atlante); k = int(np.argmin([abs(f - a.fr) for f in L["fr"]]))
    ue, ve, un, vn = L["campi"][k]; s, n_, w_, e_ = L["box"]
    verso = math.radians(a.da + 180); u0, v0 = math.sin(verso), math.cos(verso)
    spl = np.hypot(u0 * ue + v0 * un, u0 * ve + v0 * vn)
    jl = np.clip(((n_ - la_p) / (n_ - s) * L["ny"]).astype(int), 0, L["ny"] - 1)
    il = np.clip(((lo_p - w_) / (e_ - w_) * L["nx"]).astype(int), 0, L["nx"] - 1)
    l_p = spl[jl, il]
    comune = np.isfinite(c_p) & np.isfinite(l_p) & np.isfinite(r_p)
    for nome, zona in (("tutte le celle", comune), ("sottovento al Pellegrino", comune & (lo_p > 13.36) & (la_p < 38.20)),
                       ("golfo di Mondello", comune & (lo_p > 13.29) & (lo_p < 13.36) & (la_p > 38.19))):
        if zona.sum() < 10:
            print(f"{nome}: troppo poche celle ({zona.sum()})"); continue
        r = r_p[zona] / np.median(r_p[zona])
        out = [f"{nome} ({zona.sum()} celle)"]
        for etichetta, x in (("Lens", l_p[zona]), ("CFD", c_p[zona])):
            x = x / np.median(x)
            cor = np.corrcoef(x, r)[0, 1]; rms = np.sqrt(np.mean((x - r) ** 2))
            out.append(f"{etichetta}: correlazione {cor:+.2f}, scarto {rms:.3f}, escursione {np.percentile(x, 90) - np.percentile(x, 10):.2f}")
        out.append(f"radar: escursione {np.percentile(r, 90) - np.percentile(r, 10):.2f}")
        print("  ".join(out))


if __name__ == "__main__":
    main()
