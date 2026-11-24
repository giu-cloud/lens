"""Quanto vento c'è sotto costa: radar contro ECMWF e contro Lens, per distanza dalla terra, ora, forza e stagione.

  python3 costa_diagnosi.py <cartella dati> <cartella atlanti> <cartella uscita> area [area ...]
"""
import json
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_analisi as ca

BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]


def med(x):
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) >= 30 else float("nan")


def main():
    dati, atl, out = sys.argv[1:4]
    aree = sys.argv[4:]
    os.makedirs(out, exist_ok=True)
    tutto = {}
    for area in aree:
        lat, lon, scene, info = ca.carica(dati, area)
        if not scene:
            print(f"{area}: niente"); continue
        L = ca.lens(scene, lat, lon, os.path.join(atl, f"{area}.lna"), area, out)
        SU = np.array([s["su"] for s in scene]); SV = np.array([s["sv"] for s in scene])
        MU = np.array([s["mu"] for s in scene]); MV = np.array([s["mv"] for s in scene])
        sr = np.hypot(SU, SV); sm = np.hypot(MU, MV)
        sl = np.hypot(*L) if L else np.full_like(sr, np.nan)
        mai = ~np.isfinite(sr).any(axis=0)
        kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
        dist = ndimage.distance_transform_edt(~mai, sampling=(kmy, kmx))
        D = np.broadcast_to(dist, sr.shape)
        t = np.array([s["t"] for s in scene])
        ora = np.broadcast_to(((t // 3600) % 24)[:, None, None], sr.shape)
        mese = np.array([int(np.datetime64(int(x), "s").astype("datetime64[M]").astype(int) % 12 + 1) for x in t])
        estate = np.broadcast_to(((mese >= 4) & (mese <= 9))[:, None, None], sr.shape)
        ok = np.isfinite(sr) & np.isfinite(sm) & (sm > 3.0)
        rr = sr / sm; rl = sl / sm
        res = {"scene": len(scene), "celle": int(ok.sum())}
        def riga(maschera, nome):
            r = {}
            for lo, hi in BINS:
                m = ok & maschera & (D >= lo) & (D < hi)
                r[f"{lo}-{hi}"] = (med(rr[m]), med(rl[m]), int(m.sum()))
            res[nome] = r
        riga(np.ones_like(ok), "tutto")
        riga(ora < 12, "mattina (05 UTC)")
        riga(ora >= 12, "sera (17 UTC)")
        riga(estate, "aprile-settembre")
        riga(~estate, "ottobre-marzo")
        riga(sm < 6.0, "modello 3-6 m/s (6-12 nodi)")
        riga(sm >= 6.0, "modello oltre 6 m/s (12 nodi)")
        # Lens sa dove? Celle oltre 1,5 km dalla terra, divise per quanto Lens cambia il modello: che cosa dice il radar
        aff = []
        for lo, hi in [(0, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, 1.25), (1.25, 9)]:
            m = ok & (D >= 1.5) & np.isfinite(rl) & (rl >= lo) & (rl < hi)
            aff.append((f"{lo}-{hi}", med(rl[m]), med(rr[m]), int(m.sum())))
        res["affidabilita"] = aff
        # correlazione per scena tra le due forme (radar/ECMWF e Lens/ECMWF), oltre 1,5 km dalla terra
        cors = []
        for k in range(sr.shape[0]):
            m = ok[k] & (dist >= 1.5) & np.isfinite(rl[k])
            if m.sum() > 50 and np.nanstd(rl[k][m]) > 0.02:
                cors.append(np.corrcoef(rl[k][m], rr[k][m])[0, 1])
        res["correlazione_mediana"] = float(np.median(cors)) if cors else float("nan")
        res["scene_correlazione"] = len(cors)
        tutto[area] = res
        print(f"\n== {area}: {len(scene)} scene, {ok.sum()} celle con il modello oltre 3 m/s")
        print("   radar/ECMWF e Lens/ECMWF (mediane) per distanza dalla terra:")
        print("   dove Lens cambia il modello (Lens/ECMWF), che cosa dice il radar (radar/ECMWF), oltre 1,5 km dalla terra:")
        for b, a, r_, n in res["affidabilita"]:
            print(f"      Lens {b:>9}: Lens {a:.2f}  radar {r_:.2f}  ({n} celle)")
        print(f"   correlazione per scena tra le forme di Lens e del radar: mediana {res['correlazione_mediana']:.2f} su {res['scene_correlazione']} scene")
        for nome, r in res.items():
            if not isinstance(r, dict):
                continue
            print(f"   {nome:<30}" + "  ".join(f"{k:>7} km {v[0]:.2f}/{v[1]:.2f}" for k, v in r.items()))
    with open(os.path.join(out, "diagnosi.json"), "w") as f:
        json.dump(tutto, f, indent=1)


if __name__ == "__main__":
    main()
