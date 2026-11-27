"""La rotazione di Lens ha segnale? Per stazione: rotazione osservata (stazione meno modello) contro rotazione di Lens
(Lens meno modello), con modello e stazione sopra una soglia. Pendenza = quanta della rotazione di Lens è vera (1 = tutta);
errore assoluto medio e mediano del modello e di Lens; lo stesso con il peso migliore tarato fuori dalla cassaforte."""
import sys, numpy as np
from analisi import carica, dire, diff, NOMI

R = carica(sys.argv[1] if len(sys.argv) > 1 else "")
VM = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
A = np.array([r[2] for r in R]); st = np.array([r[1] for r in R]); ar = np.array([r[0] for r in R])
t, vs, ds = A[:, 0], A[:, 1], A[:, 2]
cassa = (t // (7 * 86400)) % 5 == 2
ora = ((t // 3600) % 24).astype(int)
for dove, o in (("cella della stazione", 5), ("mare aperto vicino", 15)):
    UV = {n: (A[:, o + 2 * k], A[:, o + 2 * k + 1]) for k, n in enumerate(NOMI)}
    sm = np.hypot(*UV["M"]); dM = dire(*UV["M"]); dL = dire(*UV["LA"])
    ok = (sm >= VM) & (vs >= 3.0) & np.isfinite(dL)
    oss = diff(ds, dM); rot = diff(dL, dM)
    print(f"\n== {dove}; modello >= {VM} m/s, stazione >= 3 m/s")
    print("area                 stazione            n   rot. Lens media ass.  pendenza  corr. | err. ass. medio M -> Lens | mediano M -> Lens | con il peso tarato fuori (cassaforte): M -> Lens")
    tot = []
    for k in np.unique(st[ok]):
        m = ok & (st == k)
        if m.sum() < 60:
            continue
        x = rot[m]; y = oss[m]
        pend = np.sum(x * y) / max(np.sum(x * x), 1e-9); cor = np.corrcoef(x, y)[0, 1] if np.std(x) > 0.5 else np.nan
        eM = np.abs(y); eL = np.abs(diff(y, x))
        tr = m & ~cassa; te = m & cassa
        w = np.clip(np.sum(rot[tr] * oss[tr]) / max(np.sum(rot[tr] ** 2), 1e-9), 0, 2)
        eW = np.abs(diff(oss[te], w * rot[te]))
        print(f"{ar[m][0]:20s} {k:16s} {m.sum():5d}   {np.mean(np.abs(x)):5.1f}            {pend:+5.2f}    {cor:+5.2f} |  {eM.mean():5.1f} -> {eL.mean():5.1f}   |  {np.median(eM):5.1f} -> {np.median(eL):5.1f}  | peso {w:.2f}: {np.abs(oss[te]).mean():5.1f} -> {eW.mean():5.1f} (n {te.sum()})")
        tot.append((m.sum(), eM.mean(), eL.mean(), np.abs(oss[te]).mean(), eW.mean(), te.sum(), np.mean(np.abs(x))))
    T = np.array(tot)
    print(f"media delle stazioni (ogni stazione pesa uguale): M {T[:,1].mean():.1f} -> Lens {T[:,2].mean():.1f};  cassaforte con il peso tarato: {T[:,3].mean():.1f} -> {T[:,4].mean():.1f};  migliorano {int((T[:,2] < T[:,1] - 0.5).sum())}, peggiorano {int((T[:,2] > T[:,1] + 0.5).sum())} su {len(T)}")
