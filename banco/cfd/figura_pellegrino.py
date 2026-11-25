"""Tre carte del vento attorno a Monte Pellegrino con il ponente: la carta di Pietro (Lens completo + Lens 2) a un'ora,
la CFD (OpenFOAM, aria neutra) e il radar (media dei passaggi con vento da ovest), tutte in nodi.

  python3 figura_pellegrino.py <pietro_D.json> <pietro_bin.npy> <pietro_ore.json> <ora> <caso CFD> <composito radar .npz>
        <atlante .lna> <uscita.png> [--nodi 13]

CFD e radar sono rapporti: si moltiplicano per il vento sul mare aperto della carta all'ora scelta (--nodi).
"""
import argparse
import json
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

PUNTI = {"Capo Gallo": (38.216, 13.309), "Mondello": (38.203, 13.327), "Addaura": (38.190, 13.346),
         "P. del Rotolo": (38.180, 13.372), "Arenella": (38.156, 13.380), "Acquasanta": (38.146, 13.373)}


def main():
    ap = argparse.ArgumentParser()
    for k in ("D", "bin", "ore", "ora", "caso", "radar", "atlante", "uscita"):
        ap.add_argument(k)
    ap.add_argument("--nodi", type=float, default=13.0); ap.add_argument("--da", type=float, default=278)
    ap.add_argument("--centro", default="38.170,13.353"); ap.add_argument("--zoom", default="38.12,38.235,13.29,13.44")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = json.load(open(a.D)); B = np.load(a.bin); ore = json.load(open(a.ore))["ore"]
    nx, ny = D["nx"], D["ny"]; s, n, w, e = D["interno"]
    lat = n - (np.arange(ny) + 0.5) * (n - s) / ny; lon = w + (np.arange(nx) + 0.5) * (e - w) / nx
    k = ore.index(a.ora)
    Sp = B[0, k].reshape(ny, nx) / 2.0; Dd = np.radians(B[2, k].reshape(ny, nx) * 360 / 128)
    Up, Vp = -np.sin(Dd), -np.cos(Dd)
    # CFD
    lat0, lon0 = (float(v) for v in a.centro.split(",")); G = Geometria(lat0, lon0, a.da)
    bm = open(f"{a.caso}/system/blockMeshDict").read()
    cnx, cny, cnz = (int(v) for v in re.search(r"hex \(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)", bm).groups())
    xr, yr, U10, it = cfd(a.caso, cnx, cny, cnz)
    la_c, lo_c = G.latlon(xr, yr)
    ee = U10[..., 0] * G.f[0] + U10[..., 1] * G.l[0]; nn = U10[..., 0] * G.f[1] + U10[..., 1] * G.l[1]
    vin = float(re.search(r"Uref ([\d.]+)", open(f"{a.caso}/0/U").read()).group(1))
    Sc = np.hypot(ee, nn) / vin * a.nodi
    # radar
    R = np.load(a.radar); LAr, LOr = np.meshgrid(R["lat"], R["lon"], indexing="ij"); Sr = R["rapporto"] * a.nodi
    # costa
    L = leggi_lna(a.atlante); m = L["m"]; mare = L["mare"][m:L["NY"] - m, m:L["NX"] - m]
    sb, nb, wb, eb = L["box"]
    lat_l = nb - (np.arange(L["ny"]) + 0.5) * (nb - sb) / L["ny"]; lon_l = wb + (np.arange(L["nx"]) + 0.5) * (eb - wb) / L["nx"]
    LOl, LAl = np.meshgrid(lon_l, lat_l)
    S_, N_, O_, E_ = (float(v) for v in a.zoom.split(","))
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.6), constrained_layout=True)
    LOp, LAp = np.meshgrid(lon, lat)
    pannelli = [
        (f"Carta di Pietro (Lens + Lens 2), {a.ora.replace('T', ' ')}", LOp, LAp, Sp, Up, Vp, 2),
        (f"CFD, OpenFOAM, aria neutra, vento da {a.da:.0f}°", lo_c, la_c, Sc, ee / np.maximum(np.hypot(ee, nn), 1e-6), nn / np.maximum(np.hypot(ee, nn), 1e-6), 3),
        ("Radar Sentinel-1: media di 21 passaggi con il ponente\n(vuoto entro 1,5-2 km dalla riva)", LOr, LAr, Sr, None, None, 1),
    ]
    for q, (titolo, X, Y, S, U, V, passo) in enumerate(pannelli):
        A_ = ax[q]
        im = A_.pcolormesh(X, Y, np.ma.masked_invalid(S), cmap="turbo", vmin=0, vmax=18, shading="auto")
        A_.contour(LOl, LAl, mare.astype(float), levels=[0.5], colors="black", linewidths=1.0)
        if U is not None:
            sel = (X >= O_) & (X <= E_) & (Y >= S_) & (Y <= N_)
            Xs, Ys, Us, Vs, ss = (z[::passo, ::passo] for z in (X, Y, U, V, sel))
            A_.quiver(Xs[ss], Ys[ss], Us[ss], Vs[ss], color="white", scale=48, width=0.0022, headwidth=3.5)
        for nome, (la, lo) in PUNTI.items():
            A_.plot(lo, la, "o", ms=3.5, color="black")
            A_.annotate(nome, (lo, la), xytext=(3, 3), textcoords="offset points", fontsize=8)
        A_.set_xlim(O_, E_); A_.set_ylim(S_, N_); A_.set_aspect(1 / math.cos(math.radians(0.5 * (S_ + N_))))
        A_.set_title(titolo, fontsize=10.5)
    fig.colorbar(im, ax=ax, shrink=0.85, label="nodi (CFD e radar riportati a %.0f nodi sul mare aperto)" % a.nodi)
    fig.savefig(a.uscita, dpi=120)
    print(f"-> {a.uscita} (CFD all'iterazione {it})")


if __name__ == "__main__":
    main()
