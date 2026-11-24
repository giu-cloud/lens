"""Lens tenendo solo la forma: per ogni passaggio il campo di Lens si riporta alla media del modello sul mare dell'area."""
import sys, os
import numpy as np
KN = 1.943844
AREE = ["palermo", "egadi", "eolie", "napoli", "toscano", "maddalena-bonifacio", "messina", "malta", "cagliari", "genova"]
BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]
for v in sys.argv[1:]:
    acc = {k: np.zeros((len(BINS), 3)) for k in ("mare", "lens", "lens_norm", "lens_norm_far")}
    per = {}
    for a in AREE:
        z = np.load(f"outT/{v}/campi_{a}.npz")
        R = z["radar"].astype("f8"); M = z["mare"].astype("f8"); L = z["lens_mare"].astype("f8"); dist = z["dist"]
        sea = dist > 0
        ok = np.isfinite(R) & np.isfinite(M) & np.isfinite(L) & (z["shift"] > 3.0)
        mm = np.array([np.nanmean(np.where(sea & np.isfinite(L[k]), M[k], np.nan)) for k in range(len(R))])
        ml = np.array([np.nanmean(np.where(sea & np.isfinite(L[k]), L[k], np.nan)) for k in range(len(R))])
        LN = L * (mm / ml)[:, None, None]
        far = dist > 10
        mmf = np.array([np.nanmean(np.where(far & np.isfinite(L[k]), M[k], np.nan)) if far.sum() > 30 else np.nan for k in range(len(R))])
        mlf = np.array([np.nanmean(np.where(far & np.isfinite(L[k]), L[k], np.nan)) if far.sum() > 30 else np.nan for k in range(len(R))])
        LF = L * np.where(np.isfinite(mmf / mlf), mmf / mlf, mm / ml)[:, None, None]
        D = np.broadcast_to(dist, R.shape)
        for k, F in (("mare", M), ("lens", L), ("lens_norm", LN), ("lens_norm_far", LF)):
            for i, (lo, hi) in enumerate(BINS):
                m = ok & (D >= lo) & (D < hi)
                e = (F[m] - R[m]) * KN
                acc[k][i] += [e.sum(), (e ** 2).sum(), m.sum()]
        m = ok & (D < 10)
        per[a] = {k: np.sqrt(np.mean(((F - R)[m] * KN) ** 2)) for k, F in (("mare", M), ("lens", L), ("lens_norm", LN), ("lens_norm_far", LF))}
    print(f"== {v}")
    for k, x in acc.items():
        print(f"  {k:<14}" + "".join(f"   {s / n:+5.2f} {np.sqrt(q / n):4.2f}  " for s, q, n in x))
    print("  entro 10 km, errore quadratico per area (mare / lens / lens_norm / lens_norm_far):")
    for a, r in per.items():
        print(f"    {a:<20} " + " / ".join(f"{x:.2f}" for x in r.values()))
