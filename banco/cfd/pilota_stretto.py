"""Il cancello della fase 0 (prova pilota sullo Stretto, 27/11): legge i campi della CFD stratificata (release cfd-s-stretto)
e dice, con i numeri scritti prima in lens/cfd_piano_26_11.md:
  1. convergenza: il vento a 10 m nei punti dell'asse cambia meno del 2% tra iter-300 e iter;
  2. il getto: in almeno una classe stabile, almeno 1,30 tra Messina e Reggio con il vento da nord (0°) e almeno 1,50 tra
     Messina e Villa con quello da sud (180°); la classe neutra resta intorno a 1,1 (CFD neutra del 25/11: 1,12 e 1,05);
  3. 400 m contro 150 m: la grossolana dà il getto entro 0,1 dalla fine;
  4. 15 gradi bastano: 7,5° contro la media di 0° e 15° (e 187,5° contro 180° e 195°).
Il radar (studio dello Stretto, settimane di tutti i passaggi con più di 4 m/s): 1,46-1,47 tra Messina e Reggio da nord,
1,95-1,97 tra Messina e Villa da sud.

  python3 pilota_stretto.py <cartella con i .npz> [--json uscita.json]
"""
import argparse
import glob
import json
import os
import re

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ASSE = [(38.33, 15.70), (38.29, 15.67), (38.262, 15.645), (38.245, 15.63), (38.225, 15.62), (38.20, 15.615),
        (38.175, 15.605), (38.15, 15.60), (38.12, 15.59), (38.09, 15.585), (38.06, 15.575), (38.03, 15.565), (38.00, 15.56), (37.97, 15.555)]
TRATTI = {"messina_reggio": (38.11, 38.125), "messina_villa": (38.19, 38.21), "ganzirri": (38.235, 38.255), "imbocco_nord": (38.28, 38.30)}
RADAR = {("N", "messina_reggio"): 1.47, ("S", "messina_villa"): 1.96}


def asse_fitto():
    pts = []
    for (a1, b1), (a2, b2) in zip(ASSE[:-1], ASSE[1:]):
        dd = np.hypot((a2 - a1) * 111.19, (b2 - b1) * 111.19 * np.cos(np.radians(a1)))
        for s in np.arange(0, dd, 0.25):
            pts.append((a1 + (a2 - a1) * s / dd, b1 + (b2 - b1) * s / dd))
    pts.append(ASSE[-1])
    return np.array(pts)


PTS = asse_fitto()


def profilo(p):
    z = np.load(p)
    s = np.hypot(z["u"], z["v"]); lat, lon = z["lat"], z["lon"]
    if lat[0] > lat[-1]:
        lat, s = lat[::-1], s[::-1]
    return RegularGridInterpolator((lat, lon), s, bounds_error=False)(PTS), s


def tratti(pr):
    return {k: float(np.nanmax(pr[(PTS[:, 0] >= a) & (PTS[:, 0] <= b)])) for k, (a, b) in TRATTI.items()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cartella"); ap.add_argument("--json", default="")
    a = ap.parse_args()
    casi = {}
    for p in sorted(glob.glob(os.path.join(a.cartella, "cfd_s_stretto_*m_N*_*.npz"))):
        m = re.search(r"_(\d+)m_N([\dp]+)_([\d.]+?)(_prima)?\.npz$", os.path.basename(p))
        if not m:
            continue
        dx, N, da, prima = int(m.group(1)), float(m.group(2).replace("p", ".")), float(m.group(3)), bool(m.group(4))
        casi.setdefault((dx, N, da), {})["prima" if prima else "fine"] = p
    out = {"casi": []}
    print(f"{'maglia':>6} {'N':>6} {'da':>6}  {'Messina-Reggio':>14} {'Messina-Villa':>13} {'Ganzirri':>8} {'Imbocco N':>9}  conv. (max cambio, %)")
    tab = {}
    for (dx, N, da), f in sorted(casi.items()):
        if "fine" not in f:
            print(f"{dx:6d} {N:6.3f} {da:6.1f}  manca il campo finale"); continue
        pr, s = profilo(f["fine"]); t = tratti(pr)
        conv = np.nan
        if "prima" in f:
            p0, _ = profilo(f["prima"])
            ok = np.isfinite(pr) & np.isfinite(p0)
            conv = float(np.nanmax(np.abs(pr[ok] - p0[ok]) / np.maximum(p0[ok], 0.1))) * 100
        tab[(dx, N, da)] = (pr, t)
        out["casi"].append(dict(dx=dx, N=N, da=da, **{k: round(v, 3) for k, v in t.items()}, convergenza_pct=round(conv, 2) if np.isfinite(conv) else None))
        print(f"{dx:6d} {N:6.3f} {da:6.1f}  {t['messina_reggio']:14.2f} {t['messina_villa']:13.2f} {t['ganzirri']:8.2f} {t['imbocco_nord']:9.2f}  {conv:6.2f}")
    # il cancello
    v = {}
    for dx in sorted({k[0] for k in tab}):
        stabili = [N for (d, N, da) in tab if d == dx and N > 0]
        nord = max([tab[(dx, N, 0.0)][1]["messina_reggio"] for N in set(stabili) if (dx, N, 0.0) in tab] or [np.nan])
        sud = max([tab[(dx, N, 180.0)][1]["messina_villa"] for N in set(stabili) if (dx, N, 180.0) in tab] or [np.nan])
        neu = [tab[k][1]["messina_reggio" if k[2] == 0 else "messina_villa"] for k in tab if k[0] == dx and k[1] == 0 and k[2] in (0.0, 180.0)]
        conv = [c["convergenza_pct"] for c in out["casi"] if c["dx"] == dx and c["convergenza_pct"] is not None]
        v[dx] = dict(getto_nord=nord, getto_sud=sud, neutra=neu, conv_max=max(conv) if conv else None,
                     passa_getto=bool(nord >= 1.30 and sud >= 1.50), passa_conv=bool(conv and max(conv) < 2.0))
        print(f"\nmaglia {dx} m: getto da nord {nord:.2f} (soglia 1,30, radar 1,47), da sud {sud:.2f} (soglia 1,50, radar 1,96); "
              f"neutra {', '.join(f'{x:.2f}' for x in neu)}; convergenza peggiore {v[dx]['conv_max']}% (soglia 2%)")
    if 400 in v and 150 in v:
        d = max(abs(v[400]["getto_nord"] - v[150]["getto_nord"]), abs(v[400]["getto_sud"] - v[150]["getto_sud"]))
        print(f"400 m contro 150 m: scarto massimo sul getto {d:.2f} (soglia 0,10)"); out["scarto_400_150"] = round(d, 3)
    for da, d1, d2, N in ((7.5, 0.0, 15.0, 0.01), (187.5, 180.0, 195.0, 0.01)):
        k, k1, k2 = (400, N, da), (400, N, d1), (400, N, d2)
        if k in tab and k1 in tab and k2 in tab:
            m = 0.5 * (tab[k1][0] + tab[k2][0]); e = np.abs(tab[k][0] - m)
            print(f"{da}° contro la media di {d1}° e {d2}°, lungo l'asse: scarto medio {np.nanmean(e):.3f}, massimo {np.nanmax(e):.3f} del vento al largo")
            out[f"interpolazione_{da}"] = dict(medio=round(float(np.nanmean(e)), 3), massimo=round(float(np.nanmax(e)), 3))
    out["verdetto"] = {str(k): {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in d.items()} for k, d in v.items()}
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
