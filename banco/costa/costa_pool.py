"""Tutte le aree insieme: errore medio e quadratico (nodi) per distanza dalla terra, per ciascuna lettura del modello."""
import json
import os
import sys

import numpy as np

KN = 1.943844
AREE = ["palermo", "egadi", "eolie", "napoli", "toscano", "maddalena-bonifacio", "messina", "malta", "cagliari", "genova"]
BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]


def pool(cartella, campi):
    acc = {c: [[0.0, 0.0, 0] for _ in BINS] for c in campi}
    for a in AREE:
        f = os.path.join(cartella, f"campi_{a}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        R = z["radar"].astype("f8"); D = np.broadcast_to(z["dist"], R.shape)
        ok = np.isfinite(R) & (z["shift"] > 3.0)
        for c in campi:
            ok &= np.isfinite(z[c])
        for c in campi:
            F = z[c].astype("f8")
            for i, (lo, hi) in enumerate(BINS):
                m = ok & (D >= lo) & (D < hi)
                e = (F[m] - R[m]) * KN
                acc[c][i][0] += e.sum(); acc[c][i][1] += (e ** 2).sum(); acc[c][i][2] += m.sum()
    return {c: [(s / max(n, 1), np.sqrt(q / max(n, 1)), n) for s, q, n in v] for c, v in acc.items()}


if __name__ == "__main__":
    cartella = sys.argv[1]; campi = sys.argv[2:]
    r = pool(cartella, campi)
    print(f"{'':<12}" + "".join(f"{f'{lo}-{hi} km':>16}" for lo, hi in BINS))
    for c, v in r.items():
        print(f"{c:<12}" + "".join(f"   {e:+5.2f} {q:4.2f}  " for e, q, n in v))
    print("celle", [n for _, _, n in next(iter(r.values()))])
    json.dump({c: [[float(e), float(q), int(n)] for e, q, n in v] for c, v in r.items()}, open(os.path.join(cartella, "pool.json"), "w"))
