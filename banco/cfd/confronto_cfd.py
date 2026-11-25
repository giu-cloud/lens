"""Confronto tra la CFD (caso_cfd.py, OpenFOAM) e Lens (atlante, conservazione della massa) per la stessa direzione.

  python3 confronto_cfd.py <atlante .lna> <cartella del caso> <uscita .png> --centro LAT,LON --da GRADI
        [--fr 1.0] [--zoom S,N,O,E] [--punti nome:lat:lon,...]

- CFD: vento a 10 m sul suolo, interpolato (logaritmo) tra il primo e il secondo strato di celle; riportato da x, y
  del dominio ruotato a est e nord.
- Lens: la risposta dell'atlante per la classe di stabilità più vicina a --fr, combinata per un vento unitario da --da.
- Tutte e due in rapporto al vento a 10 m sul mare all'ingresso (la CFD) o del modello (Lens).
Disegna le due carte (forza in colori, frecce della direzione) e scrive i valori nei punti scelti.
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
from caso_cfd import Geometria, leggi_atlante, leggi_lista
from controlla_atlante import leggi as leggi_lna


def leggi_campo(p, vettore):
    n, corpo = leggi_lista(p)
    if vettore:
        return np.array(corpo.replace("(", " ").replace(")", " ").split()[:3 * n], float).reshape(n, 3)
    return np.array(corpo.split()[:n], float)


def cfd(caso, nx, ny, nz):
    tempi = [d for d in os.listdir(caso) if re.fullmatch(r"\d+", d) and d != "0"]
    t = max(tempi, key=int)
    U = leggi_campo(f"{caso}/{t}/U", True).reshape(nz, ny, nx, 3)
    n, corpo = leggi_lista(f"{caso}/constant/polyMesh/points")
    P = np.array(corpo.replace("(", " ").replace(")", " ").split(), float).reshape(nz + 1, ny + 1, nx + 1, 3)
    zc = 0.25 * (P[:, :-1, :-1, 2] + P[:, 1:, :-1, 2] + P[:, :-1, 1:, 2] + P[:, 1:, 1:, 2])     # quota delle facce orizzontali
    suolo = zc[0]
    zcen = 0.5 * (zc[:-1] + zc[1:]) - suolo                                                     # quota sul suolo dei centri
    xc = 0.25 * (P[0, :-1, :-1, 0] + P[0, 1:, :-1, 0] + P[0, :-1, 1:, 0] + P[0, 1:, 1:, 0])
    yc = 0.25 * (P[0, :-1, :-1, 1] + P[0, 1:, :-1, 1] + P[0, :-1, 1:, 1] + P[0, 1:, 1:, 1])
    z1, z2 = zcen[0], zcen[1]
    w = np.clip(np.log(10.0 / z1) / np.log(z2 / z1), 0, 1)[..., None]
    U10 = U[0] * (1 - w) + U[1] * w
    return xc / 1000.0, yc / 1000.0, U10, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlante"); ap.add_argument("caso"); ap.add_argument("uscita")
    ap.add_argument("--centro", required=True); ap.add_argument("--da", type=float, required=True)
    ap.add_argument("--vento", type=float, default=8.0); ap.add_argument("--fr", type=float, default=1.0)
    ap.add_argument("--zoom", default="38.11,38.23,13.28,13.43")
    ap.add_argument("--punti", default="")
    a = ap.parse_args()
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    G = Geometria(lat0, lon0, a.da)
    bm = open(f"{a.caso}/system/blockMeshDict").read()
    nx, ny, nz = (int(v) for v in re.search(r"hex \(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)", bm).groups())
    xr, yr, U10, t = cfd(a.caso, nx, ny, nz)
    la_c, lo_c = G.latlon(xr, yr)
    e = U10[..., 0] * G.f[0] + U10[..., 1] * G.l[0]; nn = U10[..., 0] * G.f[1] + U10[..., 1] * G.l[1]
    sp_c = np.hypot(e, nn) / a.vento; dir_c = (np.degrees(np.arctan2(-e, -nn)) + 360) % 360
    # Lens: atlante
    L = leggi_lna(a.atlante)
    k = int(np.argmin([abs(f - a.fr) for f in L["fr"]]))
    ue, ve, un, vn = L["campi"][k]
    s, n_, w_, e_ = L["box"]
    lat_l = n_ - (np.arange(L["ny"]) + 0.5) * (n_ - s) / L["ny"]; lon_l = w_ + (np.arange(L["nx"]) + 0.5) * (e_ - w_) / L["nx"]
    verso = math.radians(a.da + 180); u0, v0 = math.sin(verso), math.cos(verso)
    uL = u0 * ue + v0 * un; vL = u0 * ve + v0 * vn
    sp_l = np.hypot(uL, vL); dir_l = (np.degrees(np.arctan2(-uL, -vL)) + 360) % 360
    m = L["m"]; mare = L["mare"][m:L["NY"] - m, m:L["NX"] - m]
    # valori nei punti
    righe = []
    for pt in [p for p in a.punti.split(",") if p]:
        nome, la, lo = pt.split(":"); la, lo = float(la), float(lo)
        jc = np.unravel_index(np.argmin((la_c - la) ** 2 + ((lo_c - lo) * math.cos(math.radians(la))) ** 2), la_c.shape)
        jl = int(np.argmin(abs(lat_l - la))); il = int(np.argmin(abs(lon_l - lo)))
        righe.append(dict(punto=nome, lat=la, lon=lo, cfd_forza=round(float(sp_c[jc]), 2), cfd_da=round(float(dir_c[jc])),
                          lens_forza=round(float(sp_l[jl, il]), 2), lens_da=round(float(dir_l[jl, il])), mare=bool(mare[jl, il])))
    for r in righe:
        print(f"{r['punto']:22s} CFD da {r['cfd_da']:3d}° forza {r['cfd_forza']:.2f}   Lens da {r['lens_da']:3d}° forza {r['lens_forza']:.2f}"
              f"   {'mare' if r['mare'] else 'terra'}")
    json.dump(righe, open(os.path.splitext(a.uscita)[0] + ".json", "w"), indent=1, ensure_ascii=False)
    # figura
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    S, N, O, E = (float(v) for v in a.zoom.split(","))
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.8), constrained_layout=True)
    LO, LA = np.meshgrid(lon_l, lat_l)
    for k_, (titolo, X, Y, SP, UU, VV) in enumerate((
            (f"Lens (atlante, conservazione della massa, Fr {L['fr'][k]:.1f})", LO, LA, sp_l, uL / np.maximum(sp_l, 1e-6), vL / np.maximum(sp_l, 1e-6)),
            ("CFD (OpenFOAM, RANS k-ε, aria neutra)", lo_c, la_c, sp_c, e / np.maximum(np.hypot(e, nn), 1e-6), nn / np.maximum(np.hypot(e, nn), 1e-6)))):
        A_ = ax[k_]
        im = A_.pcolormesh(X, Y, SP, cmap="viridis", vmin=0, vmax=1.3, shading="auto")
        A_.contour(LO, LA, mare.astype(float), levels=[0.5], colors="white", linewidths=1.2)
        sel = (X >= O) & (X <= E) & (Y >= S) & (Y <= N)
        passo = 3 if k_ == 0 else 3
        Xs, Ys, Us, Vs = (q[::passo, ::passo] for q in (X, Y, UU, VV)); ss = sel[::passo, ::passo]
        A_.quiver(Xs[ss], Ys[ss], Us[ss], Vs[ss], color="white", scale=45, width=0.0022, headwidth=3.5)
        for r in righe:
            A_.plot(r["lon"], r["lat"], "o", ms=4, color="orange")
            A_.annotate(r["punto"], (r["lon"], r["lat"]), xytext=(4, 4), textcoords="offset points", color="orange", fontsize=8)
        A_.set_xlim(O, E); A_.set_ylim(S, N); A_.set_aspect(1 / math.cos(math.radians(0.5 * (S + N))))
        A_.set_title(titolo, fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.8, label="forza a 10 m ÷ vento sul mare all'ingresso")
    fig.suptitle(f"Vento da {a.da:.0f}° a 10 m: Lens contro CFD", fontsize=13)
    fig.savefig(a.uscita, dpi=130)
    print(f"figura -> {a.uscita} (CFD all'iterazione {t})")


if __name__ == "__main__":
    main()
