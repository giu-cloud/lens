"""Sintesi delle varianti della lettura solo mare: errore medio e quadratico (nodi) per area e distanza, e il totale."""
import os
import sys

import numpy as np

KN = 1.943844
AREE = ["palermo", "egadi", "eolie", "napoli", "toscano", "maddalena-bonifacio", "messina", "malta", "cagliari", "genova"]
BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]
base = sys.argv[1]
varianti = sys.argv[2:]
tot = {v: {b: [0.0, 0.0, 0] for b in BINS} for v in ["shift"] + varianti}
print(f"{'':<22}" + "".join(f"{f'{lo}-{hi} km':>16}" for lo, hi in BINS))
for a in AREE:
    righe = {}
    for v in varianti:
        f = os.path.join(base, v, f"campi_{a}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        R = z["radar"].astype("f8"); S0 = z["shift"].astype("f8"); M = z["mare"].astype("f8"); D = np.broadcast_to(z["dist"], R.shape)
        ok = np.isfinite(R) & np.isfinite(S0) & np.isfinite(M) & (S0 > 3.0)
        for nome, F in ([("shift", S0)] if "shift" not in righe else []) + [(v, M)]:
            r = []
            for b in BINS:
                m = ok & (D >= b[0]) & (D < b[1])
                if m.sum() < 50:
                    r.append(None); continue
                e = F[m] - R[m]
                r.append((np.mean(e) * KN, np.sqrt(np.mean(e ** 2)) * KN))
                t = tot[nome][b]; t[0] += e.sum() * KN; t[1] += (e ** 2).sum() * KN * KN; t[2] += m.sum()
            righe[nome] = r
    print(f"\n{a}")
    for nome, r in righe.items():
        print(f"  {nome:<20}" + "".join(f"{'':>16}" if x is None else f"   {x[0]:+5.1f} {x[1]:4.1f}   " for x in r))
print("\nTUTTE LE AREE")
for nome, t in tot.items():
    print(f"  {nome:<20}" + "".join(f"   {t[b][0] / max(t[b][2], 1):+5.2f} {np.sqrt(t[b][1] / max(t[b][2], 1)):4.2f}  " for b in BINS))
