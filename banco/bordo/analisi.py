"""Il bordo dei riquadri di Lens (28/11): la prova sui riquadri finti.

Per ogni area e lato, il lato del riquadro è spostato verso l'interno di 15 km; il livello del radar è rifatto senza
niente oltre quel lato. Nella striscia tra il lato finto e quello vero si sa che cosa c'è (il radar), e si confrontano:
  M         il modello da solo;
  oggi      Lens dentro, il modello fuori, raccordo di 1,5 miglia verso l'interno (LensWind della 0.59);
  alone H   dentro uguale, ma il raccordo va verso il modello corretto dal rapporto di Lens sul bordo, e fuori il modello
            resta corretto da quel rapporto che si spegne in H km (lisciato lungo il bordo su 5 km);
  C pieno   Lens con il riquadro intero (quello che si avrebbe con un riquadro più grande: il meglio possibile).
  python3 analisi.py [aree]
"""
import sys, numpy as np

L = "/home/claude/lavoro/bordo"; KN = 1 / 0.5144; B = 1.5 * 1.852; RI = 15.0
BOX = {"palermo": (37.945, 38.287, 13.108, 13.612), "messina": (37.90, 38.35, 15.40, 15.75), "eolie": (38.00, 38.86, 14.25, 15.30),
       "monaco": (43.45, 43.85, 7.15, 8.00), "corsica-nord": (42.45, 43.10, 8.25, 9.65), "dalmazia": (42.85, 43.55, 15.95, 17.25),
       "cagliari": (38.85, 39.25, 8.85, 9.65)}
AREE = sys.argv[1:] or list(BOX)
FASCE = [(-15, -10), (-10, -5), (-5, -2), (-2, 0), (0, B), (B, 8)]
HS = (5.0, 10.0, 15.0)


def lato(X, box, s):
    S, N, W, E = box; ky = 111.19; kx = 111.19 * np.cos(np.radians((S + N) / 2))
    la, lo = X["lat"], X["lon"]
    dn, ds, dw, de = (N - la) * ky, (la - S) * ky, (lo - W) * kx, (E - lo) * kx      # km dai lati veri
    if s == "n": d, x, altri = dn - RI, lo * kx, np.minimum.reduce([ds, dw, de])
    elif s == "s": d, x, altri = ds - RI, lo * kx, np.minimum.reduce([dn, dw, de])
    elif s == "o": d, x, altri = dw - RI, la * ky, np.minimum.reduce([dn, ds, de])
    else: d, x, altri = de - RI, la * ky, np.minimum.reduce([dn, ds, dw])
    return d, x, altri


def previsioni(X, s, d, x):
    mod = X["mod"]; r = X["r_" + s]; rf = X["r_full"]
    w = np.clip(d / B, 0, 1)
    P = {"M": mod.copy(), "C pieno": mod * rf}
    P["oggi"] = np.where(d <= 0, mod, mod * (w * r + (1 - w)))
    # il rapporto sul bordo, per passaggio, lisciato lungo il bordo (gaussiana di 5 km) dalla striscia interna di 3 km
    R = np.ones(len(mod))
    for t in np.unique(X["t"]):
        k = X["t"] == t; st = k & (d > 0) & (d < 3) & np.isfinite(r)
        if st.sum() < 3:
            continue
        xs, ms, rs = x[st], mod[st], r[st]
        idx = np.where(k)[0]
        wg = np.exp(-0.5 * ((x[idx][:, None] - xs[None, :]) / 5.0) ** 2)
        num = (wg * (ms * rs)[None, :]).sum(1); den = (wg * ms[None, :]).sum(1)
        R[idx] = np.where(den > 1e-6, np.clip(num / np.maximum(den, 1e-6), 0.5, 2.0), 1.0)
    for H in HS:
        P[f"alone {H:.0f}"] = np.where(d <= 0, mod * (1 + (R - 1) * np.clip(1 + d / H, 0, 1)), mod * (w * r + (1 - w) * R))
    return P, R


def main():
    global righe
    righe = []
    for a in AREE:
        try:
            Z = np.load(f"{L}/bordo_{a}.npz"); X = {k: Z[k] for k in Z.files if k != "box"}
        except FileNotFoundError:
            continue
        ok0 = np.isfinite(X["vero"]) & (X["vero"] > 0.5) & (X["mod"] > 2.0)
        for s in "nsoe":
            d, x, altri = lato(X, BOX[a], s)
            k = ok0 & (d > -RI) & (d < 8) & (altri > 3) & np.isfinite(X["r_" + s]) & np.isfinite(X["r_full"])
            if k.sum() < 200:
                continue
            P, R = previsioni({kk: v[k] for kk, v in X.items()}, s, d[k], x[k])
            righe.append(dict(area=a, lato=s, d=d[k], t=X["t"][k], vero=X["vero"][k], P=P, R=R, dist=X["dist"][k]))

    nomi = list(righe[0]["P"])
    print(f"{len(righe)} lati in {len(set(r['area'] for r in righe))} aree; {sum(len(r['d']) for r in righe)} righe")
    print("\nErrore medio / quadratico (nodi) per fascia di distanza dal lato finto (negativo = fuori, dove oggi c'è il modello):")
    print("fascia km     n     " + "  ".join(f"{n:>13s}" for n in nomi))
    for lo_, hi in FASCE:
        e = {n: np.concatenate([(r["P"][n] - r["vero"])[(r["d"] >= lo_) & (r["d"] < hi)] for r in righe]) * KN for n in nomi}
        print(f"{lo_:5.1f}..{hi:4.1f} {len(e['M']):6d}  " + "  ".join(f"{e[n].mean():+5.2f}/{np.sqrt((e[n] ** 2).mean()):4.2f}" for n in nomi))

    # per area: fuori (da -10 a 0) e nel raccordo (da 0 a 2,8): quadratico, e differenza alone 10 - oggi con intervallo
    print("\nPer area, da 10 km fuori al lato finto: quadratico di M, oggi, alone 10, C pieno; alone 10 meno oggi (intervallo 90%, ricampionando i passaggi)")
    rng = np.random.default_rng(4)
    for a in list(dict.fromkeys(r["area"] for r in righe)) + ["TUTTE"]:
        rr = [r for r in righe if a == "TUTTE" or r["area"] == a]
        sel = [(r, (r["d"] >= -10) & (r["d"] < 0)) for r in rr]
        e = {n: np.concatenate([((r["P"][n] - r["vero"]) * KN)[m] for r, m in sel]) for n in ("M", "oggi", "alone 10", "C pieno")}
        tt = np.concatenate([r["t"][m] + (hash(r["area"] + r["lato"]) % 997) * 1e10 for r, m in sel])
        U = np.unique(tt); bb = []
        pos = {u: np.where(tt == u)[0] for u in U}
        for _ in range(300):
            p = np.concatenate([pos[u] for u in rng.choice(U, len(U))])
            bb.append(np.sqrt((e["alone 10"][p] ** 2).mean()) - np.sqrt((e["oggi"][p] ** 2).mean()))
        q = lambda n: np.sqrt((e[n] ** 2).mean())   # noqa: E731
        print(f"  {a:13s} n {len(e['M']):6d}  {q('M'):.2f}  {q('oggi'):.2f}  {q('alone 10'):.2f}  {q('C pieno'):.2f}   {q('alone 10') - q('oggi'):+.3f} ({np.percentile(bb, 5):+.3f} … {np.percentile(bb, 95):+.3f})")

    # il gradiente finto: per lato e passaggio, il profilo medio lungo il bordo a passi di 1 km; la pendenza più forte tra -6 e +4 km
    print("\nIl gradiente sul bordo: pendenza più forte del profilo medio lungo il bordo, tra 6 km fuori e 4 km dentro (nodi per miglio), mediana e 90° percentile")
    bins = np.arange(-6, 5, 1.0)
    for n in ("oggi", "alone 5", "alone 10", "alone 15", "C pieno"):
        sl = []
        for r in righe:
            for t in np.unique(r["t"]):
                k = r["t"] == t
                prof = [np.mean(r["P"][n][k & (r["d"] >= b0) & (r["d"] < b0 + 1)]) if np.any(k & (r["d"] >= b0) & (r["d"] < b0 + 1)) else np.nan for b0 in bins]
                prof = np.array(prof) * KN
                dd = np.abs(np.diff(prof)) * 1.852
                if np.isfinite(dd).sum() >= 6:
                    sl.append(np.nanmax(dd))
        sl = np.array(sl)
        print(f"  {n:9s}: {np.median(sl):.2f}  {np.percentile(sl, 90):.2f}   (su {len(sl)} profili)")


if __name__ == "__main__":
    main()
