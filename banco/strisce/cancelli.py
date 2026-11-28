"""I cancelli della prova pilota delle strisce (lens/strisce_piano_28_11.md), dalle finestre di analizza.py.

  python3 cancelli.py <area> [soglia di qualità]

Cancello 1: mare aperto (oltre 20 km): quota di finestre buone e differenza mediana dal modello; qui si fissa la soglia.
Cancello 2: le stazioni: radar contro stazione, modello contro stazione, alla stessa ora, finestre entro 6 km.
Cancello 3: ripetibilità: la mappa della rotazione media (radar meno ECMWF) in due metà dei passaggi.
"""
import json, os, sys
import numpy as np

L = "/home/claude/lavoro"; D = "/home/claude/dati"
area = sys.argv[1]
Z = np.load(f"{L}/strisce/finestre_{area}.npz"); area = area.replace("-buoni", ""); X = {k: Z[k] for k in Z.files}
dif = lambda a, b: (a - b + 180) % 360 - 180
passaggi = np.unique(X["t"])
print(f"{area}: {len(passaggi)} passaggi, {len(X['t'])} finestre")

# ---------- cancello 1
lar = X["dist"] > 20
print("\nCANCELLO 1, mare aperto (oltre 20 km dalla costa): soglia di qualità, quota di finestre tenute, |radar - ECMWF| mediano e medio")
scelta = None
for thr in (0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50):
    m = lar & (X["q"] >= thr)
    if m.sum() < 50:
        continue
    med = np.median(np.abs(X["dif"][m])); quota = m.sum() / lar.sum()
    print(f"  q >= {thr:.2f}: {quota * 100:5.1f}%   mediana {med:5.1f}   media {np.mean(np.abs(X['dif'][m])):5.1f}")
    if scelta is None and med <= 15.0:
        scelta = thr
if len(sys.argv) > 2:
    scelta = float(sys.argv[2])
if scelta is None:
    scelta = 0.25
m = lar & (X["q"] >= scelta); quota = m.sum() / lar.sum(); med = np.median(np.abs(X["dif"][m]))
print(f"  soglia fissata: q >= {scelta:.2f} (la più bassa con mediana entro 15 gradi). Quota {quota * 100:.0f}%, mediana {med:.1f} gradi.")
# per passaggio
pp = []
for t in passaggi:
    k = m & (X["t"] == t); a = lar & (X["t"] == t)
    if a.sum() >= 50:
        pp.append((k.sum() / a.sum(), np.median(np.abs(X["dif"][k])) if k.sum() >= 20 else np.nan, np.median(X["ecm_v"][a])))
pp = np.array(pp)
print(f"  per passaggio: quota di finestre buone mediana {np.median(pp[:, 0]) * 100:.0f}% (da {pp[:, 0].min() * 100:.0f}% a {pp[:, 0].max() * 100:.0f}%); passaggi con mediana entro 20 gradi: {int(np.sum(pp[:, 1] <= 20))} su {int(np.sum(np.isfinite(pp[:, 1])))}")
# le cuciture: quante direzioni buone cadono entro 5 gradi dagli assi dell'immagine (attese a caso: 11%)
a_img = X["ang_img"][m]; assi = (np.minimum(a_img % 90, 90 - a_img % 90) < 5).mean()
print(f"  direzioni entro 5 gradi dagli assi dell'immagine: {assi * 100:.0f}% (a caso sarebbero l'11%)")
ic = m & np.isfinite(X["icon"])
if ic.sum() > 100:
    print(f"  dove c'è ICON-2I ({ic.sum()} finestre): |radar - ICON-2I| mediana {np.median(np.abs(dif(X['sar'][ic], X['icon'][ic]))):.1f}; |radar - ECMWF| {np.median(np.abs(X['dif'][ic])):.1f}; |ICON-2I - ECMWF| {np.median(np.abs(dif(X['icon'][ic], X['ecm'][ic]))):.1f}")
ok1 = quota >= 0.30 and med <= 20.0
print("  ESITO:", "passa" if ok1 else "NON passa")

# ---------- cancello 2
E = [s for s in json.load(open(f"{L}/stazioni_elenco.json")) if s["area"] == area and not (s.get("quota") or 0) > 200]
print("\nCANCELLO 2, le stazioni (finestre di mare entro 6 km, qualità sopra la soglia, stazione ad almeno 3 m/s)")
righe = []
for s in E:
    p = next((f"{D}/{d}/{s['fonte']}_{s['codice']}.npz" for d in ("stazioni_r4", "stazioni_r3") if os.path.exists(f"{D}/{d}/{s['fonte']}_{s['codice']}.npz")), None)
    if p is None:
        continue
    S = np.load(p); o = np.argsort(S["t"]); ts, vs, ds = S["t"][o], S["vento"][o], S["direzione"][o]
    km = 111.19 * np.hypot(X["lat"] - s["lat"], (X["lon"] - s["lon"]) * np.cos(np.radians(s["lat"])))
    for t in passaggi:
        k = (X["t"] == t) & (km <= 6) & (X["q"] >= scelta)
        j = np.searchsorted(ts, t); best = None
        for q in (j - 1, j):
            if 0 <= q < len(ts) and abs(ts[q] - t) <= 1200 and np.isfinite(vs[q]) and np.isfinite(ds[q]) and (best is None or abs(ts[q] - t) < abs(ts[best] - t)):
                best = q
        if best is None or vs[best] < 3.0 or k.sum() < 3:
            continue
        cm = lambda a: np.degrees(np.angle(np.mean(np.exp(1j * np.radians(a))))) % 360
        kk = (X["t"] == t) & (km <= 6)
        righe.append((s["nome"][:14], t, ds[best], cm(X["sar"][k]), cm(X["ecm"][kk]), cm(X["icon"][kk]) if np.isfinite(X["icon"][kk]).all() else np.nan, k.sum()))
if righe:
    nomi = np.array([r[0] for r in righe]); R = np.array([r[1:] for r in righe], float)
    eS = dif(R[:, 2], R[:, 1]); eM = dif(R[:, 3], R[:, 1]); eI = dif(R[:, 4], R[:, 1])
    print("  stazione          n | errore ass. medio: radar  ECMWF  ICON-2I | medio con segno: radar  ECMWF")
    for nm in list(np.unique(nomi)) + ["TUTTE"]:
        k = (nomi == nm) if nm != "TUTTE" else np.ones(len(nomi), bool)
        ki = k & np.isfinite(eI)
        print(f"  {nm:16s} {k.sum():3d} | {np.abs(eS[k]).mean():5.1f}  {np.abs(eM[k]).mean():5.1f}  {np.abs(eI[ki]).mean() if ki.sum() else np.nan:5.1f} ({ki.sum()}) | {eS[k].mean():+6.1f} {eM[k].mean():+6.1f}")
    g = np.abs(eM).mean() - np.abs(eS).mean()
    rng = np.random.default_rng(5); tt = np.unique(R[:, 0]); bb = []
    for _ in range(500):
        pick = np.concatenate([np.where(R[:, 0] == x)[0] for x in rng.choice(tt, len(tt))])
        bb.append(np.abs(eM[pick]).mean() - np.abs(eS[pick]).mean())
    print(f"  guadagno del radar sul modello: {g:+.1f} gradi (intervallo 90%: {np.percentile(bb, 5):+.1f} … {np.percentile(bb, 95):+.1f}), su {len(R)} confronti")
    kc = np.array(["Corse" in x or "CORSE" in x.upper() for x in nomi])
    if kc.sum() >= 5:
        print(f"  Cap Corse: rotazione media rispetto a ECMWF: stazione {(-eM[kc]).mean():+.1f}, radar {dif(R[kc, 2], R[kc, 3]).mean():+.1f} (n {kc.sum()})")
    ok2 = len(R) >= 20 and g >= 5.0
    print("  ESITO:", "passa" if ok2 else "NON passa")
else:
    print("  nessun confronto")

# ---------- cancello 3
print("\nCANCELLO 3, ripetibilità della rotazione (radar meno ECMWF), celle di 0,02 gradi entro 15 km dalla costa")
m = (X["q"] >= scelta) & (X["dist"] <= 15)
cel = (np.floor(X["lat"] / 0.02) * 10000 + np.floor(X["lon"] / 0.02)).astype(np.int64)
# il settore di vento più frequente (ECMWF medio del passaggio, in quattro settori)
sett = {t: int(((np.degrees(np.angle(np.mean(np.exp(1j * np.radians(X["ecm"][X["t"] == t]))))) % 360) + 45) // 90) % 4 for t in passaggi}
conta = np.bincount(list(sett.values()), minlength=4); nomi_s = ["nord", "est", "sud", "ovest"]
print("  passaggi per settore del vento:", {nomi_s[k]: int(conta[k]) for k in range(4)})
for nome_s, quali in [("tutti i venti", list(passaggi))] + [(nomi_s[k], [t for t in passaggi if sett[t] == k]) for k in np.argsort(-conta)[:2]]:
    if len(quali) < 8:
        continue
    rng = np.random.default_rng(3); cc = []
    for _ in range(200):
        p = rng.permutation(quali); A = set(p[: len(p) // 2]); B = set(p[len(p) // 2:])
        def media(S_):
            k = m & np.isin(X["t"], list(S_)); ids, inv = np.unique(cel[k], return_inverse=True)
            n = np.bincount(inv); z = np.bincount(inv, np.cos(np.radians(2 * X["dif"][k]))) + 1j * np.bincount(inv, np.sin(np.radians(2 * X["dif"][k])))
            return dict(zip(ids[n >= 6], (0.5 * np.degrees(np.angle(z)))[n >= 6]))
        a, b = media(A), media(B); com = sorted(set(a) & set(b))
        if len(com) >= 30:
            cc.append(np.corrcoef([a[x] for x in com], [b[x] for x in com])[0, 1])
    if cc:
        print(f"  {nome_s:14s} ({len(quali)} passaggi): correlazione tra le due metà, mediana {np.median(cc):+.2f} (da {np.percentile(cc, 10):+.2f} a {np.percentile(cc, 90):+.2f})")
