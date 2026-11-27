"""Quanto dell'errore di direzione è "fisso" (dipende solo dalla direzione del vento al largo, cioè da rilievo e canali) e
quanto no: per stazione, la correzione media per settore di 30 gradi imparata fuori dalla cassaforte e provata dentro.
È il tetto di qualunque modello del rilievo fermo come Lens. In più: per ora del giorno e per forza."""
import sys, numpy as np
from analisi import carica, dire, diff, NOMI
R = carica(""); VM = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
A = np.array([r[2] for r in R]); st = np.array([r[1] for r in R]); ar = np.array([r[0] for r in R])
t, vs, ds = A[:, 0], A[:, 1], A[:, 2]; cassa = (t // (7 * 86400)) % 5 == 2; ora = ((t // 3600) % 24).astype(int)
o = 15
UV = {n: (A[:, o + 2 * k], A[:, o + 2 * k + 1]) for k, n in enumerate(NOMI)}
sm = np.hypot(*UV["M"]); dM = dire(*UV["M"]); dL = dire(*UV["LA"])
ok = (sm >= VM) & (vs >= 3.0) & np.isfinite(dL)
oss = diff(ds, dM); rot = diff(dL, dM); sett = (dM // 30).astype(int)
print("stazione: n in cassaforte | errore ass. medio: modello, Lens, modello + correzione fissa per settore, Lens + correzione fissa del suo residuo")
T = []
for k in np.unique(st[ok]):
    m = ok & (st == k); tr = m & ~cassa; te = m & cassa
    if te.sum() < 25 or tr.sum() < 80:
        continue
    def corr(res):
        c = np.zeros(12)
        for s in range(12):
            q = tr & (sett == s)
            if q.sum() >= 5:
                c[s] = np.degrees(np.arctan2(np.mean(np.sin(np.radians(res[q]))), np.mean(np.cos(np.radians(res[q])))))
        return c
    cM = corr(oss); resL = diff(oss, rot); cL = corr(resL)
    eM = np.abs(oss[te]).mean(); eL = np.abs(resL[te]).mean()
    eMc = np.abs(diff(oss[te], cM[sett[te]])).mean(); eLc = np.abs(diff(resL[te], cL[sett[te]])).mean()
    print(f"{ar[m][0]:20s} {k:16s} {te.sum():4d} | {eM:5.1f}  {eL:5.1f}  {eMc:5.1f}  {eLc:5.1f}")
    T.append((eM, eL, eMc, eLc))
T = np.array(T); print("media delle stazioni:           | " + "  ".join(f"{x:5.1f}" for x in T.mean(0)))
print("\nper ora UTC (tutte le stazioni, modello >= %.1f): errore ass. medio del modello / di Lens / n" % VM)
for h in range(0, 24, 3):
    m = ok & (ora == h); print(f"  {h:02d}: {np.abs(oss[m]).mean():5.1f} / {np.abs(diff(oss[m], rot[m])).mean():5.1f} / {m.sum()}")
print("per forza del modello (m/s), tutte le righe con stazione >= 2 m/s:")
ok2 = (vs >= 2.0) & np.isfinite(dL)
for a_, b_ in ((2.5, 4), (4, 6), (6, 8), (8, 11), (11, 99)):
    m = ok2 & (sm >= a_) & (sm < b_); print(f"  {a_}-{b_}: {np.abs(oss[m]).mean():5.1f} / {np.abs(diff(oss[m], rot[m])).mean():5.1f} / {m.sum()}")
