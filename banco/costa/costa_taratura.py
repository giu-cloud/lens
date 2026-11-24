"""Lens contro il radar con pesi diversi dei suoi pezzi (ridosso, fascia costiera, rilievo), sopra ECMWF letto dai soli
punti di mare. Per ogni variante e area: errore medio e quadratico (nodi) per distanza dalla terra, e la correlazione
per scena tra la forma di Lens e quella del radar (rapporto col modello), oltre 1,5 km.

  python3 costa_taratura.py <uscita> variante [variante ...]     (varianti in VARIANTI)
"""
import json
import os
import subprocess
import sys

import numpy as np

S = os.path.dirname(os.path.abspath(__file__))
KN = 1.943844
AREE = ["palermo", "egadi", "eolie", "napoli", "toscano", "maddalena-bonifacio", "messina", "malta", "cagliari", "genova"]
VARIANTI = {
    "lens": {},
    "senza_ridosso": {"shift.lens.wake": 0},
    "senza_costa": {"shift.lens.coast": 0},
    "solo_rilievo": {"shift.lens.wake": 0, "shift.lens.coast": 0},
    "meta": {"shift.lens.wake": 0.3, "shift.lens.coast": 0.5},
    "stabile_corto": {"shift.lens.sfacMax": 1.0},
    "niente": {"shift.lens.wake": 0, "shift.lens.coast": 0, "shift.lens.terrain": 0},
    "rilievo_meta": {"shift.lens.wake": 0, "shift.lens.coast": 0, "shift.lens.terrain": 0.5},
    "ridosso_terzo": {"shift.lens.wake": 0.2, "shift.lens.coast": 0, "shift.lens.sfacMax": 1.0},
    "costa_meta": {"shift.lens.coast": 0.5},
    "radar": {"shift.lens.wake": 0.3, "shift.lens.sfacMax": 1.0, "shift.lens.coast": 0, "shift.lens.terrain": 0.5},
    "solo_ridosso_lieve": {"shift.lens.wake": 0.2, "shift.lens.sfacMax": 1.0, "shift.lens.coast": 0, "shift.lens.terrain": 0},
    "radar_lieve": {"shift.lens.wake": 0.2, "shift.lens.sfacMax": 1.0, "shift.lens.coast": 0, "shift.lens.terrain": 0.25},
}
BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]


def main():
    out = sys.argv[1]
    for v in sys.argv[2:]:
        d = os.path.join(out, v)
        os.makedirs(d, exist_ok=True)
        env = dict(os.environ, LENS_PROPS=json.dumps(VARIANTI[v]), UT=f"{S}/../ut2")
        subprocess.run(["python3", f"{S}/costa_ecmwf.py", f"{S}/dati", f"{S}/atlanti", f"{S}/ecmwf", d, *AREE], env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print(f"{v}: fatto", flush=True)


def riassunto(out, varianti):
    tab = {}
    for v in varianti:
        for a in AREE:
            f = os.path.join(out, v, f"campi_{a}.npz")
            if not os.path.exists(f):
                continue
            z = np.load(f)
            R = z["radar"].astype("f8"); M = z["mare"].astype("f8"); L = z["lens_mare"].astype("f8"); dist = z["dist"]
            ok = np.isfinite(R) & np.isfinite(M) & np.isfinite(L) & (M > 3.0)
            D = np.broadcast_to(dist, R.shape)
            r = {}
            for lo, hi in BINS:
                m = ok & (D >= lo) & (D < hi)
                if m.sum() < 50:
                    continue
                r[f"{lo}-{hi}"] = dict(err=float(np.mean(L[m] - R[m]) * KN), rmse=float(np.sqrt(np.mean((L[m] - R[m]) ** 2)) * KN),
                                       err_mare=float(np.mean(M[m] - R[m]) * KN), rmse_mare=float(np.sqrt(np.mean((M[m] - R[m]) ** 2)) * KN), n=int(m.sum()))
            cs = []
            for k in range(R.shape[0]):
                m = ok[k] & (dist >= 1.5)
                if m.sum() > 50:
                    a_ = L[k][m] / M[k][m]; b_ = R[k][m] / M[k][m]
                    if np.std(a_) > 0.02:
                        cs.append(np.corrcoef(a_, b_)[0, 1])
            r["corr"] = float(np.median(cs)) if cs else float("nan")
            tab.setdefault(a, {})[v] = r
    return tab


if __name__ == "__main__":
    if sys.argv[1] == "--riassunto":
        out = sys.argv[2]; vs = sys.argv[3:]
        tab = riassunto(out, vs)
        json.dump(tab, open(os.path.join(out, "taratura.json"), "w"), indent=1)
        for a, rv in tab.items():
            print(f"\n== {a}   (errore medio / quadratico in nodi; 'mare' = ECMWF solo mare senza Lens)")
            any_v = next(iter(rv.values()))
            print(f"   {'':<16}" + "".join(f"{b:>14}" for b in any_v if b != "corr") + "   corr")
            print(f"   {'mare':<16}" + "".join(f"  {x['err_mare']:+5.1f} {x['rmse_mare']:4.1f}  " for b, x in any_v.items() if b != "corr"))
            for v, r in rv.items():
                print(f"   {v:<16}" + "".join(f"  {x['err']:+5.1f} {x['rmse']:4.1f}  " for b, x in r.items() if b != "corr") + f"  {r['corr']:+.2f}")
    else:
        main()
