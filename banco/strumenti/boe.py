"""Prova di Lens 2 contro le boe, a tutte le ore: ECMWF e ICON-2I, grezzi e corretti con le librerie del radar.

  python3 boe.py <boa .npz> <area> <librerie> <ecmwf passi 3-6> <ecmwf passi 9-12> <stab 3-6> <925 9-12> <icon2i ore> <rilievo> <uscita.json>

- ECMWF: corse delle 00 e delle 12, passi 3, 6, 9 e 12 ore (valide alle 03, 06, 09, 12 e 15, 18, 21, 00 UTC), letto nelle
  celle della libreria dai soli punti di mare a 0,25° (un nodo è di terra se il rilievo nel punto è sopra zero, come la
  costa di Shift). Il vento a 925 hPa è la media dei nodi dell'area.
- ICON-2I: storia di Open-Meteo, tutte le ore, letto dai soli punti di mare della sua griglia.
- La boa: misure entro 15 minuti dall'ora del modello. Il vento della boa è a pochi metri sul mare: per non confondere
  la quota con la correzione, ogni previsione si giudica anche dopo averla moltiplicata per un solo numero (lo stesso
  per tutte le ore), stimato su tutta la serie.
"""
import datetime as dt
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from banco_icon2i import terra_nodi
from libreria_radar import Libreria
from rilievo import Rilievo
from stabilita import leggi

KN = 1.943844


def lettura_ecmwf(lat, lon, terra, lat_c, lon_c):
    """Pesi bilineari dai soli nodi di mare per le celle (vettori lat_c, lon_c)."""
    dla = lat[1] - lat[0]; dlo = lon[1] - lon[0]
    out = []
    for la, lo in zip(lat_c, lon_c):
        fj = (la - lat[0]) / dla; j0 = int(np.floor(fj)); y = fj - j0
        fi = (lo - lon[0]) / dlo; i0 = int(np.floor(fi)); x = fi - i0
        pts = []
        for j, i, w in ((j0, i0, (1 - x) * (1 - y)), (j0, i0 + 1, x * (1 - y)), (j0 + 1, i0, (1 - x) * y), (j0 + 1, i0 + 1, x * y)):
            if 0 <= j < len(lat) and 0 <= i < len(lon) and not terra[j, i] and w > 0:
                pts.append((j, i, w))
        if not pts:
            best = None; bd = 1e9
            for j in range(j0 - 2, j0 + 4):
                for i in range(i0 - 2, i0 + 4):
                    if 0 <= j < len(lat) and 0 <= i < len(lon) and not terra[j, i]:
                        dd = (j - fj) ** 2 + (i - fi) ** 2
                        if dd < bd:
                            bd, best = dd, (j, i, 1.0)
            pts = [best] if best else []
        s = sum(w for _, _, w in pts)
        out.append([(j, i, w / s) for j, i, w in pts])
    return out


def main():
    boa_p, area, cart_l, e36, e912, s36, p912, icon_p, ril_p, usc = sys.argv[1:11]
    B = np.load(boa_p)
    ril = Rilievo(ril_p)
    LE = Libreria(os.path.join(cart_l, f"lra1_ecmwf_{area}.lra"))
    LI_p = os.path.join(cart_l, f"lra1_icon-2i_{area}.lra")
    LI = Libreria(LI_p) if os.path.exists(LI_p) else None
    lat_c = np.array([LE.lat0 + (c // LE.nx) * LE.dlat for c in LE.celle])
    lon_c = np.array([LE.lon0 + (c % LE.nx) * LE.dlon for c in LE.celle])
    blat, blon = float(np.nanmedian(B["lat"])), float(np.nanmedian(B["lon"]))
    cb = int(np.argmin((lat_c - blat) ** 2 + ((lon_c - blon) * np.cos(np.radians(blat))) ** 2))
    print(f"boa {blat:.3f} {blon:.3f}: cella della libreria a {111 * np.sqrt((lat_c[cb] - blat) ** 2 + ((lon_c[cb] - blon) * np.cos(np.radians(blat))) ** 2):.2f} km, "
          f"a {LE.dist[cb]:.1f} km da terra", flush=True)
    # ICON-2I
    if LI is not None and os.path.exists(icon_p):
        Z = np.load(icon_p); ilat, ilon = Z["lat"], Z["lon"]; it = Z["t"]
        o = np.argsort(it); it = it[o]; IU = Z["u"][o]; IV = Z["v"][o]
        it_terra = terra_nodi(ril, ilat, ilon)
        lat_ci = np.array([LI.lat0 + (c // LI.nx) * LI.dlat for c in LI.celle]); lon_ci = np.array([LI.lon0 + (c % LI.nx) * LI.dlon for c in LI.celle])
        pesi_i = lettura_ecmwf(ilat, ilon, it_terra, lat_ci, lon_ci)
        cbi = int(np.argmin((lat_ci - blat) ** 2 + ((lon_ci - blon) * np.cos(np.radians(blat))) ** 2))
    else:
        Z = None
    righe = []
    pesi = None; nodi_area = None
    bt = B["t"]; bv = B["vento"]
    for k, t in enumerate(bt):
        if not np.isfinite(bv[k]):
            continue
        tv = int(round(t / 3600.0)) * 3600
        if abs(t - tv) > 900:
            continue
        d = dt.datetime.fromtimestamp(tv, dt.timezone.utc)
        if d.hour % 3:
            continue
        run = d.replace(hour=0 if 3 <= d.hour <= 12 else 12, minute=0, second=0)
        if d.hour == 0:
            run = run - dt.timedelta(days=1)
        passo = int((d - run).total_seconds() // 3600)
        nome = f"ecmwf_{run:%Y%m%d%H}.grb2"
        f_w = os.path.join(e36, nome) if passo <= 6 else os.path.join(e912, nome.replace("ecmwf_", "e912_"))
        f_s = os.path.join(s36, nome.replace("ecmwf_", "stab_")) if passo <= 6 else os.path.join(p912, nome.replace("ecmwf_", "p912_"))
        if not os.path.exists(f_w):
            continue
        lat, lon, C = leggi(f_w)
        if (("10u", passo)) not in C:
            continue
        if pesi is None:
            LA, LO = np.meshgrid(lat, lon, indexing="ij")
            terra = np.zeros(LA.shape, bool)
            for j in range(LA.shape[0]):
                for i in range(LA.shape[1]):
                    if abs(LA[j, i] - blat) < 1.5 and abs(LO[j, i] - blon) < 2:
                        _, _, hh = ril.griglia(0, LA[j, i] - 0.003, LA[j, i] + 0.003, LO[j, i] - 0.004, LO[j, i] + 0.004)
                        terra[j, i] = hh.max() > 0
            pesi = lettura_ecmwf(lat, lon, terra, lat_c, lon_c)
            s_, n_, w_, e_ = lat_c.min(), lat_c.max(), lon_c.min(), lon_c.max()
            nodi_area = (LA >= s_ - 0.25) & (LA <= n_ + 0.25) & (LO >= w_ - 0.25) & (LO <= e_ + 0.25)
        U = C[("10u", passo)]; V = C[("10v", passo)]
        u = np.array([sum(w * U[j, i] for j, i, w in p) for p in pesi], "f4"); v = np.array([sum(w * V[j, i] for j, i, w in p) for p in pesi], "f4")
        a925 = b925 = None
        if os.path.exists(f_s):
            _, _, S = leggi(f_s)
            if ("u@925", passo) in S:
                u9 = float(S[("u@925", passo)][nodi_area].mean()); v9 = float(S[("v@925", passo)][nodi_area].mean())
                s9 = np.hypot(u9, v9); d9 = np.arctan2(-u9, -v9); a925, b925 = s9 * np.sin(d9), s9 * np.cos(d9)
        largo = LE.dist >= 8
        if largo.sum() < 20:
            largo = LE.dist >= np.quantile(LE.dist[LE.dist > 0], 0.75)
        S0 = float(np.median(np.hypot(u, v)[largo])); D0 = np.arctan2(-u[largo].mean(), -v[largo].mean())
        a0, b0 = S0 * np.sin(D0), S0 * np.cos(D0)
        if a925 is None:
            a925, b925 = a0, b0
        sera = 1.0 if 12 <= d.hour <= 23 else 0.0
        dirc = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
        G, SP = LE.fattore(a0, b0, a925, b925, sera, dirc, dispersione=True)
        r = {"t": int(tv), "ora": d.hour, "boa": float(bv[k]), "ecmwf": float(np.hypot(u[cb], v[cb])), "ecmwf_c": float(np.hypot(u[cb], v[cb]) * G[cb]),
             "disp": float(SP[cb]), "dir": float(dirc[cb]), "passo": passo}
        if Z is not None:
            kk = np.searchsorted(it, tv)
            if kk < len(it) and it[kk] == tv:
                Uk = IU[kk].astype("f4"); Vk = IV[kk].astype("f4")
                ui = np.array([sum(w * Uk[j, i] for j, i, w in pp) for pp in pesi_i], "f4")
                vi = np.array([sum(w * Vk[j, i] for j, i, w in pp) for pp in pesi_i], "f4")
                largo_i = LI.dist >= 8
                if largo_i.sum() < 20:
                    largo_i = LI.dist >= np.quantile(LI.dist[LI.dist > 0], 0.75)
                S0i = float(np.median(np.hypot(ui, vi)[largo_i])); D0i = np.arctan2(-ui[largo_i].mean(), -vi[largo_i].mean())
                dici = (np.degrees(np.arctan2(-ui, -vi)) + 360) % 360
                Gi = LI.fattore(S0i * np.sin(D0i), S0i * np.cos(D0i), a925, b925, sera, dici)
                r["icon"] = float(np.hypot(ui[cbi], vi[cbi])); r["icon_c"] = float(r["icon"] * Gi[cbi])
        righe.append(r)
    json.dump(righe, open(usc, "w"))
    print(f"{len(righe)} confronti -> {usc}")


if __name__ == "__main__":
    main()
