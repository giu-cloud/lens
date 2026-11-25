"""Vento a 925 hPa della scena da una corsa con un anticipo dato (per le prove a uno o due giorni).

  python3 vento925.py <cartella p*_*.grb2> <prefisso dei file> <anticipo ore> <cartella del metro> <rilievo> <uscita>

Come stabilita.py: media sui nodi di mare entro 0,25° dal riquadro e di terra entro 0,5°; la corsa è quella con
l'anticipo più vicino a quello chiesto (entro 6 ore), interpolata tra i suoi due passi. Scrive v925_<area>.npz con
t, v925, d925 (provenienza).
"""
import datetime as dt
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rilievo import Rilievo
from stabilita import AREE, leggi


def main():
    cart, pre, ant, metro, ril_p, usc = sys.argv[1:7]
    ant = float(ant)
    os.makedirs(usc, exist_ok=True)
    ril = Rilievo(ril_p)
    disp = set(os.listdir(cart))
    cache = {}
    for area in AREE:
        p = os.path.join(metro, f"campi_{area}.npz")
        if not os.path.exists(p):
            continue
        z = np.load(p); t = z["t"]
        la0, la1, lo0, lo1 = z["lat"].min(), z["lat"].max(), z["lon"].min(), z["lon"].max()
        v9 = np.full(len(t), np.nan, "f4"); d9 = np.full(len(t), np.nan, "f4")
        mask = None
        for s, ts in enumerate(t):
            d = dt.datetime.fromtimestamp(int(ts), dt.timezone.utc).replace(tzinfo=None)
            scelta = None
            for back in range(0, 10):
                base = d.replace(minute=0, second=0, microsecond=0)
                r = base.replace(hour=(base.hour // 12) * 12) - dt.timedelta(hours=12 * back)
                lead = (d - r).total_seconds() / 3600
                if lead < ant - 6:
                    continue
                if lead > ant + 6:
                    break
                n = f"{pre}_{r:%Y%m%d%H}.grb2"
                if n in disp:
                    scelta = (n, lead); break
            if scelta is None:
                continue
            n, lead = scelta
            if n not in cache:
                if len(cache) > 40:
                    cache.clear()
                cache[n] = leggi(os.path.join(cart, n))
            lat, lon, C = cache[n]
            passi = sorted({st for (k, st) in C if k == "u@925"})
            if not passi:
                continue
            if mask is None:
                LA, LO = np.meshgrid(lat, lon, indexing="ij")
                terra = np.zeros(LA.shape, bool)
                for j in range(LA.shape[0]):
                    for i in range(LA.shape[1]):
                        if la0 - 0.6 <= LA[j, i] <= la1 + 0.6 and lo0 - 0.6 <= LO[j, i] <= lo1 + 0.6:
                            _, _, hh = ril.griglia(0, LA[j, i] - 0.04, LA[j, i] + 0.04, LO[j, i] - 0.05, LO[j, i] + 0.05)
                            terra[j, i] = (hh > 0).mean() > 0.5
                vicino = (LA >= la0 - 0.25) & (LA <= la1 + 0.25) & (LO >= lo0 - 0.25) & (LO <= lo1 + 0.25)
                largo = (LA >= la0 - 0.5) & (LA <= la1 + 0.5) & (LO >= lo0 - 0.5) & (LO <= lo1 + 0.5)
                mask = (vicino & ~terra) | (largo & terra)
            a, b = passi[0], passi[-1]
            w = 0.0 if b == a else min(max((lead - a) / (b - a), 0.0), 1.0)
            u = (1 - w) * C[("u@925", a)] + w * C[("u@925", b)]; v = (1 - w) * C[("v@925", a)] + w * C[("v@925", b)]
            um, vm = float(u[mask].mean()), float(v[mask].mean())
            v9[s] = np.hypot(um, vm); d9[s] = (np.degrees(np.arctan2(-um, -vm)) + 360) % 360
        np.savez(os.path.join(usc, f"v925_{area}.npz"), t=t, v925=v9, d925=d9)
        print(f"{area}: {np.isfinite(v9).sum()} su {len(t)}", flush=True)


if __name__ == "__main__":
    main()
