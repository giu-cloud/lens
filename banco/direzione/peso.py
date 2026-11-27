"""Peso e tetto della rotazione di Lens sul mare vicino alla stazione: tarati fuori dalla cassaforte, provati dentro."""
import sys, numpy as np
from analisi import carica, dire, diff, NOMI
R = carica(""); VM = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
A = np.array([r[2] for r in R]); st = np.array([r[1] for r in R])
t, vs, ds = A[:, 0], A[:, 1], A[:, 2]; cassa = (t // (7 * 86400)) % 5 == 2
for dove, o in (("cella", 5), ("mare", 15)):
    UV = {n: (A[:, o + 2 * k], A[:, o + 2 * k + 1]) for k, n in enumerate(NOMI)}
    sm = np.hypot(*UV["M"]); dM = dire(*UV["M"]); dL = dire(*UV["LA"])
    ok = (sm >= VM) & (vs >= 3.0) & np.isfinite(dL)
    oss = diff(ds, dM); rot = diff(dL, dM)
    S = [k for k in np.unique(st[ok]) if (ok & (st == k)).sum() >= 60]
    def punteggio(w, c, mk):
        return np.mean([np.mean(np.abs(diff(oss[mk & (st == k)], np.clip(w * rot[mk & (st == k)], -c, c)))) for k in S if (mk & (st == k)).sum() >= 10])
    print(f"\n{dove}: {len(S)} stazioni; errore assoluto medio (media delle stazioni), fuori dalla cassaforte / cassaforte")
    print("  modello:            %.2f / %.2f" % (punteggio(0, 1, ok & ~cassa), punteggio(0, 1, ok & cassa)))
    for w in (0.5, 0.75, 1.0, 1.25, 1.5):
        print(f"  peso {w:.2f}: " + "  ".join(f"tetto {c:3d}: {punteggio(w, c, ok & ~cassa):.2f}/{punteggio(w, c, ok & cassa):.2f}" for c in (10, 20, 30, 180)))
