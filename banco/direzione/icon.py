"""Direzione alle stazioni: ECMWF (come lo riceve Lens), Lens su ECMWF, ICON-2I (2 km, punto più vicino e media a 3 km),
sulle stesse ore. E il tetto della correzione fissa per settore sopra ICON-2I."""
import json, sys, numpy as np
from analisi import carica, dire, diff, NOMI
R = carica(""); VM = 5.0
E = {s["fonte"] + "_" + s["codice"]: s for s in json.load(open("/home/claude/lavoro/stazioni_elenco.json"))}
A = np.array([r[2] for r in R]); st = np.array([r[1] for r in R]); ar = np.array([r[0] for r in R])
t, vs, ds = A[:, 0], A[:, 1], A[:, 2]; cassa = (t // (7 * 86400)) % 5 == 2
o = 15; UV = {n: (A[:, o + 2 * k], A[:, o + 2 * k + 1]) for k, n in enumerate(NOMI)}
sm = np.hypot(*UV["M"]); dM = dire(*UV["M"]); dL = dire(*UV["LA"])
dI = np.full(len(t), np.nan); sI = np.full(len(t), np.nan)
for a in np.unique(ar):
    try: Z = np.load(f"/home/claude/dati/icon2i_3h/icon2i_{a}.npz")
    except FileNotFoundError: continue
    it = Z["t"]; lat = Z["lat"]; lon = Z["lon"]; U = Z["u"]; V = Z["v"]
    for k in np.unique(st[ar == a]):
        s = E[k]; j = int(np.argmin(np.abs(lat - s["lat"]))); i = int(np.argmin(np.abs(lon - s["lon"])))
        if abs(lat[j] - s["lat"]) > 0.05 or abs(lon[i] - s["lon"]) > 0.05: continue
        m = np.where((ar == a) & (st == k))[0]; q = np.searchsorted(it, t[m]); q = np.clip(q, 0, len(it) - 1); buono = it[q] == t[m]
        u = U[q, j, i].astype(float); v = V[q, j, i].astype(float)
        dI[m[buono]] = dire(u, v)[buono]; sI[m[buono]] = np.hypot(u, v)[buono]
ok = (sm >= VM) & (vs >= 3) & np.isfinite(dI) & np.isfinite(dL)
print("stazione            n | err. ass. medio: ECMWF  Lens su ECMWF  ICON-2I | medio con segno: ECMWF ICON-2I | cassaforte, con la correzione fissa per settore: ECMWF ICON-2I")
T = []
for k in np.unique(st[ok]):
    m = ok & (st == k)
    if m.sum() < 60: continue
    eM = diff(ds[m], dM[m]); eL = diff(ds[m], dL[m]); eI = diff(ds[m], dI[m])
    def tetto(dmod):
        tr = m & ~cassa; te = m & cassa; res = diff(ds, dmod); sett = (dmod // 30).astype(int); c = np.zeros(12)
        for s_ in range(12):
            qq = tr & (sett == s_)
            if qq.sum() >= 5: c[s_] = np.degrees(np.arctan2(np.mean(np.sin(np.radians(res[qq]))), np.mean(np.cos(np.radians(res[qq])))))
        return np.abs(diff(res[te], c[sett[te]])).mean() if te.sum() >= 15 else np.nan
    row = (np.abs(eM).mean(), np.abs(eL).mean(), np.abs(eI).mean(), eM.mean(), eI.mean(), tetto(dM), tetto(np.nan_to_num(dI)))
    print(f"{ar[m][0][:12]:12s} {k:16s} {m.sum():4d} | {row[0]:5.1f}  {row[1]:5.1f}  {row[2]:5.1f} | {row[3]:+5.1f} {row[4]:+5.1f} | {row[5]:5.1f} {row[6]:5.1f}")
    T.append(row)
T = np.array(T); print("media delle stazioni:               | " + "  ".join(f"{x:5.1f}" for x in np.nanmean(T, 0)))

# la direzione dell'insieme: media dei vettori ECMWF e ICON-2I con peso w per ICON-2I; e la rotazione di Lens sopra
uM, vM = UV["M"]; rot = diff(dL, dM)
S = [k for k in np.unique(st[ok]) if (ok & (st == k)).sum() >= 60]
def media(e, mk):
    return np.mean([np.abs(e[mk & (st == k)]).mean() for k in S if (mk & (st == k)).sum() >= 10])
print("\ndirezione dell'insieme con peso w a ICON-2I (vettori di forza unitaria): fuori dalla cassaforte / cassaforte; poi lo stesso con la rotazione di Lens (calcolata su ECMWF) aggiunta, intera e a metà")
for w in (0.0, 0.4, 0.6, 0.8, 1.0):
    a = np.radians(dM); b = np.radians(np.nan_to_num(dI))
    d = np.degrees(np.arctan2((1 - w) * np.sin(a) + w * np.sin(b), (1 - w) * np.cos(a) + w * np.cos(b))) % 360
    e0 = diff(ds, d); e1 = diff(ds, d + rot); e5 = diff(ds, d + 0.5 * rot)
    print(f"  w {w:.1f}: {media(e0, ok & ~cassa):.1f} / {media(e0, ok & cassa):.1f}   con Lens: {media(e1, ok & ~cassa):.1f} / {media(e1, ok & cassa):.1f}   con mezza: {media(e5, ok & ~cassa):.1f} / {media(e5, ok & cassa):.1f}")
