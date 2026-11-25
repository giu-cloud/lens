"""La CFD di un sito contro il radar, passaggio per passaggio: quanto vale il disegno della CFD rispetto a Lens 0.54.

  python3 taratura_radar.py --sito palermo --area palermo --cfd <cartella con cfd_<sito>_XXX.X.npz>
        --mare <atlante .lna o terreno .npz> --campi <campi_<area>.npz> --tabella <tabella6.npz>
        --vero <radar_vero_tabella6.npy> --uscita <risultati .json> [--npz <righe .npz>]

Per ogni cella del radar (circa 1,3 km) la forza della CFD è la media delle sue celle di mare (150 m) dentro la cella,
per ognuna delle 16 direzioni. In ogni riga della tabella (una cella in un passaggio) la direzione è quella del modello
al largo (D0) e la CFD si interpola tra le due direzioni vicine. Il rapporto della CFD è la forza nella cella divisa per
la mediana sulle celle al largo (più di 8 km da terra), così si confronta con il modello letto dal mare in quella cella.

Bersaglio: radar (vento vero) diviso il modello letto dal mare nella cella. Forme provate, in logaritmo:
  lens      c + b log(Lens 0.54)
  cfd       c + a log(CFD)
  insieme   c + b log(Lens 0.54) + a log(CFD)
Le costanti si stimano fuori dalla cassaforte e si provano sulle settimane della cassaforte ((t // 7 giorni) % 5 == 2).
Anche per terzi di forza del vento al largo e di stabilità (dT850), per vedere se la forza della CFD cambia.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

AREE14 = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova",
          "pontine", "corsica-nord", "monaco", "cicladi"]
SETTIMANA = 7 * 86400


def pila_cfd(cartella, sito):
    """Le direzioni calcolate (ordinate) e i campi u, v (k, ny, nx), con lat e lon della griglia."""
    voci = []
    for p in glob.glob(os.path.join(cartella, f"cfd_{sito}_*.npz")):
        if re.search(r"_v[\d.]+\.npz$", p):
            continue
        z = np.load(p)
        voci.append((float(z["da"]), z["u"], z["v"], z["lat"], z["lon"]))
    voci.sort(key=lambda x: x[0])
    da = np.array([v[0] for v in voci])
    U = np.stack([v[1] for v in voci]); V = np.stack([v[2] for v in voci])
    return da, U, V, voci[0][3], voci[0][4]


def mare_su(fonte, lat, lon):
    """Maschera del mare sulla griglia (lat, lon) dal terreno del sito o dall'atlante (il più vicino)."""
    if fonte.endswith(".npz"):
        z = np.load(fonte); la, lo, m = z["lat"], z["lon"], z["mare"]
    else:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from caso_cfd import leggi_atlante
        A = leggi_atlante(fonte); la, lo, m = A["lat"], A["lon"], A["mare"]
    j = np.clip(np.round((la[0] - lat) / (la[0] - la[1])).astype(int), 0, len(la) - 1)
    i = np.clip(np.round((lon - lo[0]) / (lo[1] - lo[0])).astype(int), 0, len(lo) - 1)
    return m[np.ix_(j, i)]


def su_celle_radar(S, mare, lat, lon, rlat, rlon):
    """Media delle celle di mare valide della CFD dentro ogni cella del radar: (k, nr, nc) e frazione coperta."""
    k = S.shape[0]
    dla = abs(rlat[1] - rlat[0]); dlo = abs(rlon[1] - rlon[0])
    out = np.full((k, len(rlat), len(rlon)), np.nan); cop = np.zeros((len(rlat), len(rlon)))
    for jr, a in enumerate(rlat):
        rj = np.where(np.abs(lat - a) <= dla / 2)[0]
        if not len(rj):
            continue
        for ir, b in enumerate(rlon):
            ri = np.where(np.abs(lon - b) <= dlo / 2)[0]
            if not len(ri):
                continue
            blocco = S[:, rj[0]:rj[-1] + 1, ri[0]:ri[-1] + 1]
            mm = mare[rj[0]:rj[-1] + 1, ri[0]:ri[-1] + 1]
            ok = np.isfinite(blocco).all(axis=0) & mm
            cop[jr, ir] = ok.mean()
            if ok.mean() >= 0.5:
                out[:, jr, ir] = blocco[:, ok].mean(axis=1)
    return out, cop


def interpola(da, F, theta):
    """Le due direzioni calcolate attorno a theta (gradi da cui viene il vento) e il peso della seconda."""
    from stazioni_cfd import vicini
    return vicini(da, theta)


def adatta(X, y):
    X1 = np.column_stack([np.ones(len(y))] + X)
    c, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return c


def prevedi(c, X):
    return c[0] + sum(ci * xi for ci, xi in zip(c[1:], X))


def misure(obs, prev, dist, scena, cella):
    """Errore quadratico, errore sistematico (per cella: media dell'errore su tutti i passaggi, poi la radice della media
    dei quadrati sulle celle, come nello studio statistico), gli stessi entro 3 km, e la correlazione del disegno
    (mediana sui passaggi della correlazione tra le celle)."""
    e = prev - obs
    def sist(mask):
        if not mask.any():
            return None
        ids, inv = np.unique(cella[mask], return_inverse=True)
        me = np.bincount(inv, e[mask]) / np.bincount(inv)
        n = np.bincount(inv)
        return float(np.sqrt(np.sum(n * me ** 2) / n.sum()))
    vic = dist < 3
    corr = []
    for s in np.unique(scena):
        m = scena == s
        if m.sum() >= 15 and np.std(obs[m]) > 0 and np.std(prev[m]) > 0:
            corr.append(np.corrcoef(obs[m], prev[m])[0, 1])
    tutto = np.ones(len(e), bool)
    return dict(n=int(len(e)), rmse=float(np.sqrt(np.mean(e ** 2))), sist=sist(tutto), media=float(np.mean(e)),
                rmse_3km=float(np.sqrt(np.mean(e[vic] ** 2))) if vic.any() else None, sist_3km=sist(vic), n_3km=int(vic.sum()),
                corr_disegno=float(np.median(corr)) if corr else None, passaggi=len(corr))


def main():
    ap = argparse.ArgumentParser()
    for k in ("sito", "area", "cfd", "mare", "campi", "tabella", "vero", "uscita"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--npz", default="")
    a = ap.parse_args()
    da, U, V, lat, lon = pila_cfd(a.cfd, a.sito)
    print(f"{a.sito}: {len(da)} direzioni, griglia {U.shape[1]} x {U.shape[2]}")
    S = np.hypot(U, V)
    mare = mare_su(a.mare, lat, lon)
    C = np.load(a.campi); rlat, rlon, dist_c = C["lat"].astype(float), C["lon"].astype(float), C["dist"]
    Sr, cop = su_celle_radar(S, mare, lat, lon, rlat, rlon)
    Soff = []
    for k in range(len(da)):
        v = np.nan
        for soglia in (8, 6, 4):          # nello Stretto le celle a più di 8 km da terra sono poche
            sel = (dist_c > soglia) & np.isfinite(Sr[k])
            if sel.sum() >= 5:
                v = np.nanmedian(Sr[k][sel]); break
        Soff.append(v)
    Soff = np.array(Soff)
    print("celle del radar coperte dalla CFD:", int(np.isfinite(Sr).all(axis=0).sum()), "di", Sr[0].size,
          "· mediana al largo per direzione:", np.round(Soff, 3).tolist())
    T = np.load(a.tabella)
    ia = AREE14.index(a.area)
    m = T["area"] == ia
    j, i = T["j"][m].astype(int), T["i"][m].astype(int)
    vero = np.load(a.vero)[m]
    mod = T["mod"][m]; D0 = T["D0"][m]; S0 = T["S0"][m]; dist = T["dist"][m]; t = T["t"][m]; scena = T["scena"][m]
    qL = T["q_lens054"][m]; dT850 = T["dT850"][m]; dT925 = T["dT925"][m]
    k0, k1, w = interpola(da, Sr, D0)
    s = (1 - w) * Sr[k0, j, i] + w * Sr[k1, j, i]
    so = (1 - w) * Soff[k0] + w * Soff[k1]
    rc = s / so
    con_lens = np.isfinite(qL).mean() > 0.5          # per le aree nuove Lens 0.54 non c'è nella tabella
    if not con_lens:
        qL = np.ones_like(qL); print("Lens 0.54 non c'è nella tabella per quest'area: solo modello e CFD")
    ok = np.isfinite(rc) & np.isfinite(vero) & (vero > 0.5) & (mod > 2.0) & np.isfinite(qL) & (qL > 0) & (S0 > 2.0)
    print(f"righe: {m.sum()}, valide con la CFD: {ok.sum()}")
    cass = ((t // SETTIMANA) % 5 == 2)
    y = np.log(vero / mod)
    lL, lC = np.log(np.where(qL > 0, qL, np.nan)), np.log(np.where(rc > 0, rc, np.nan))
    forme = {"modello": [], "lens": [lL], "cfd": [lC], "insieme": [lL, lC]} if con_lens else {"modello": [], "cfd": [lC]}
    tr, te = ok & ~cass, ok & cass
    ris = {"sito": a.sito, "area": a.area, "direzioni": da.tolist(), "righe_addestramento": int(tr.sum()), "righe_cassaforte": int(te.sum()),
           "celle_coperte": int(np.isfinite(Sr).all(axis=0).sum()), "forme": {}}
    for nome, X in forme.items():
        c = adatta([x[tr] for x in X], y[tr]) if X else np.array([np.mean(y[tr])])
        p = prevedi(c, [x[te] for x in X]) if X else np.full(te.sum(), c[0])
        obs_v = vero[te]; prev_v = mod[te] * np.exp(p)
        cel = (j[te] * 1000 + i[te])
        mis = misure(obs_v, prev_v, dist[te], scena[te], cel)
        mis_r = misure(y[te], p, dist[te], scena[te], cel)
        ris["forme"][nome] = dict(coefficienti=[float(x) for x in c], nodi={k: (v / 0.5144 if isinstance(v, float) and k.startswith(("rmse", "sist", "media")) else v) for k, v in mis.items()},
                                  corr_disegno_log=mis_r["corr_disegno"])
        k3 = lambda v: v / 0.5144 if v is not None else float("nan")
        print(f"{nome:8s} coeff {np.round(c, 3).tolist()} · cassaforte: rmse {k3(mis['rmse']):.2f} kn, sist {k3(mis['sist']):.2f}, "
              f"entro 3 km rmse {k3(mis['rmse_3km']):.2f} sist {k3(mis['sist_3km']):.2f} (n {mis['n_3km']}) · disegno {mis_r['corr_disegno'] if mis_r['corr_disegno'] is None else round(mis_r['corr_disegno'], 3)}")
    # la forza della CFD per terzi di vento al largo e di stabilità (sull'insieme di tutte le righe valide, forma "insieme")
    strati = {}
    for nome, q in (("vento_largo", S0), ("dT850", dT850), ("dT925", dT925)):
        lim = np.nanpercentile(q[ok], [33.3, 66.7])
        for b, (lo_, hi_) in enumerate(((-np.inf, lim[0]), (lim[0], lim[1]), (lim[1], np.inf))):
            mm = ok & (q > lo_) & (q <= hi_)
            if mm.sum() < 200:
                continue
            c = adatta([lL[mm], lC[mm]], y[mm]) if con_lens else np.array([np.nan, np.nan, np.nan])
            c1 = adatta([lC[mm]], y[mm])
            strati[f"{nome}_{b}"] = dict(intervallo=[float(lo_), float(hi_)], n=int(mm.sum()), insieme=[float(x) for x in c], solo_cfd=[float(x) for x in c1])
            print(f"  {nome} {b} ({lo_:.1f}, {hi_:.1f}] n {mm.sum()}: insieme b {c[1]:.2f} a {c[2]:.2f} · solo CFD a {c1[1]:.2f}")
    ris["strati"] = strati
    json.dump(ris, open(a.uscita, "w"), indent=1)
    if a.npz:
        np.savez_compressed(a.npz, ok=ok, cass=cass, y=y, lL=lL, lC=lC, mod=mod, vero=vero, dist=dist, scena=scena, D0=D0, S0=S0,
                            j=j, i=i, w=w, dT850=dT850, Sr=Sr, Soff=Soff, da=da, rlat=rlat, rlon=rlon)


if __name__ == "__main__":
    main()
