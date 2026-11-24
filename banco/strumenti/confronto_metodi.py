"""Confronto di metodi per correggere il vento del modello vicino alla costa, provati su aree che non hanno visto.

  python3 confronto_metodi.py <tabella.npz> <uscita> [metodi separati da virgola] [--bersaglio rel|ass] [--cassaforte]

Il metodo stima un fattore g per cella e passaggio: vento corretto = g × modello letto dal mare (× k0 nel bersaglio
relativo, dove k0 è il rapporto radar/modello al largo della stessa scena: così si guarda solo la forma vicino alla
costa, non la taratura generale del modello contro il radar).
Prova: si lascia fuori un'area alla volta (si allena sulle altre nove, si giudica su quella). Le settimane della
cassaforte (una su cinque) non si usano mai, se non con --cassaforte, alla fine.
Misure sull'area lasciata fuori:
  - errore medio e quadratico in nodi, per distanza dalla terra;
  - errore sistematico: per ogni cella e ottante di direzione con almeno 5 passaggi, lo scarto medio; se ne fa lo
    scarto quadratico tolto il rumore del radar (la varianza dentro il gruppo divisa per il numero di passaggi);
  - scarti medi per classi del rilievo (ridosso, valle, canale, terra davanti).
"""
import json
import os
import sys
import time

import numpy as np

KN = 1.943844
AREE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]
DBIN = [0, 1.5, 3, 5, 10, 99]
SETTIMANA = 7 * 86400


def carica(p):
    T = dict(np.load(p))
    T["y"] = (T["radar"] / T["mod"]).astype("f4")
    for k in [k for k in T if k.startswith("q_")]:
        T[k] = np.nan_to_num(T[k], nan=1.0)
    T["lmin"] = np.minimum(T["lsx"], T["ldx"])
    T["canale"] = np.minimum(T["lsx"] + T["ldx"], 60.0)
    T["sera"] = (T["ora"] >= 12).astype("f4")
    T["msin"] = np.sin(2 * np.pi * T["mese"] / 12).astype("f4")
    T["mcos"] = np.cos(2 * np.pi * T["mese"] / 12).astype("f4")
    T["rS0"] = (T["mod"] / np.maximum(T["S0"], 0.5)).astype("f4")
    T["ott"] = (np.round(T["dir"] / 45.0).astype(int) % 8).astype("i1")
    T["cassa"] = ((T["t"] // SETTIMANA) % 5 == 2)
    return T


STATICHE = ["dist", "terra", "h15", "ang", "hc", "valle", "f5", "f15", "sotto", "hsotto", "angsotto", "lmin", "canale",
            "angw", "fw15", "h10", "hm10", "terra10"]
SCENA = ["S0", "rS0", "dDir", "sera", "msin", "mcos"]


# ---------------------------------------------------------------- metodi: fit(T, idx, bersaglio, peso) -> predici(T, idx)

def m_modello(T, tr, y, w):
    return lambda te: np.ones(len(te), dtype="f4")


def m_lens054(T, tr, y, w):
    return lambda te: np.nan_to_num(T["q_lens054"][te], nan=1.0)


def m_costante(T, tr, y, w):
    c = np.sum(w * y) / np.sum(w)
    return lambda te: np.full(len(te), c, dtype="f4")


def _binned(chiavi_tr, chiavi_te, y, w, minimo=2000):
    k_all = np.concatenate([chiavi_tr, chiavi_te])
    u, inv = np.unique(k_all, return_inverse=True)
    itr, ite = inv[:len(chiavi_tr)], inv[len(chiavi_tr):]
    sw = np.bincount(itr, weights=w, minlength=len(u)); swy = np.bincount(itr, weights=w * y, minlength=len(u))
    cnt = np.bincount(itr, minlength=len(u))
    glob = np.sum(w * y) / np.sum(w)
    val = np.where(cnt >= minimo, swy / np.maximum(sw, 1e-9), glob)
    return val[ite].astype("f4")


def m_distanza(T, tr, y, w):
    b = np.array([0, 0.75, 1.5, 2.25, 3, 4, 5, 6.5, 8, 10, 12.5, 15, 20, 99])
    def chi(ix): return np.digitize(T["dist"][ix], b)
    ktr = chi(tr)
    return lambda te: _binned(ktr, chi(te), y, w)


def m_dist_rilievo(T, tr, y, w):
    """distanza × tipo di vento: dal mare, da terra bassa, da terra alta (angolo del rilievo)."""
    b = np.array([0, 1.5, 3, 5, 8, 12, 99])
    def chi(ix):
        d = np.digitize(T["dist"][ix], b)
        tipo = np.where(T["terra"][ix] >= 29.9, 0, np.where(T["ang"][ix] < 1, 1, np.where(T["ang"][ix] < 4, 2, 3)))
        return d * 10 + tipo
    ktr = chi(tr)
    return lambda te: _binned(ktr, chi(te), y, w)


def _basi(T, ix):
    d = T["dist"][ix]; terra = T["terra"][ix]; da_terra = (terra < 29.9).astype("f4")
    cols = [np.exp(-d / 1.5), np.exp(-d / 5),
            da_terra * np.exp(-terra / 3), da_terra * np.exp(-terra / 10),
            np.minimum(T["ang"][ix], 15) / 10, np.minimum(T["angw"][ix], 20) / 10,
            np.log1p(T["h15"][ix]) / 7 * da_terra, np.log1p(T["hc"][ix]) / 6 * da_terra,
            np.clip(T["valle"][ix], -300, 1000) / 500, T["f5"][ix], T["f15"][ix], T["fw15"][ix],
            np.exp(-T["sotto"][ix] / 3), np.minimum(T["angsotto"][ix], 15) / 10,
            np.exp(-T["lmin"][ix] / 3) * (T["canale"][ix] < 20), np.exp(-T["canale"][ix] / 10),
            T["terra10"][ix], T["h10"][ix] / 1000, T["hm10"][ix] / 300]
    return np.stack(cols, axis=1).astype("f8")


def m_lineare(T, tr, y, w):
    X = _basi(T, tr)
    Xs = X * np.sqrt(w)[:, None]
    A = Xs.T @ Xs + 1e-3 * np.sum(w) * np.eye(X.shape[1]) / len(w)
    c = np.sum(w * y) / np.sum(w)
    b = np.linalg.solve(A, Xs.T @ ((y - c) * np.sqrt(w)))
    return lambda te: (c + _basi(T, te) @ b).astype("f4")


def _hgb(caratt, monot=None, **kw):
    def fit(T, tr, y, w):
        from sklearn.ensemble import HistGradientBoostingRegressor
        rng = np.random.default_rng(0)
        sub = tr if len(tr) <= 1_200_000 else np.sort(rng.choice(tr, 1_200_000, replace=False))
        pos = np.searchsorted(tr, sub)
        X = np.stack([T[c][sub] for c in caratt], axis=1)
        par = dict(max_iter=250, learning_rate=0.06, max_leaf_nodes=31, min_samples_leaf=3000, l2_regularization=1.0,
                   early_stopping=False, random_state=0)
        par.update(kw)
        if monot:
            par["monotonic_cst"] = [monot.get(c, 0) for c in caratt]
        g = HistGradientBoostingRegressor(**par)
        g.fit(X, y[pos], sample_weight=w[pos])
        return lambda te: g.predict(np.stack([T[c][te] for c in caratt], axis=1)).astype("f4")
    return fit


METODI = {
    "modello": m_modello,
    "lens054": m_lens054,
    "costante": m_costante,
    "distanza": m_distanza,
    "dist_rilievo": m_dist_rilievo,
    "lineare": m_lineare,
    "alberi_rilievo": _hgb(STATICHE),
    "alberi_tutto": _hgb(STATICHE + SCENA),
    "alberi_tutto_lens": _hgb(STATICHE + SCENA + ["q_lens054"]),
}


def aggiungi_fisica(T):
    nomi = [k for k in T if k.startswith("q_") and k != "q_lens054" and np.isfinite(T[k]).all()]
    if nomi:
        METODI["alberi_tutto_fisica"] = _hgb(STATICHE + SCENA + nomi)
        METODI["alberi_rilievo_fisica"] = _hgb(STATICHE + nomi)
    return nomi


# ---------------------------------------------------------------- misure

def misure(T, ix, g, bers):
    R = T["radar"][ix]; base = T["mod"][ix] * (T["k0"][ix] if bers == "rel" else 1.0)
    P = g * base
    e = (P - R) * KN
    d = T["dist"][ix]
    out = {"n": int(len(ix)), "media": float(np.mean(e)), "rmse": float(np.sqrt(np.mean(e * e)))}
    for lo, hi in zip(DBIN[:-1], DBIN[1:]):
        m = (d >= lo) & (d < hi)
        if m.sum() > 100:
            out[f"d{lo}-{hi}"] = [float(np.mean(e[m])), float(np.sqrt(np.mean(e[m] ** 2))), int(m.sum())]
    # sistematico per cella e ottante
    key = (T["area"][ix].astype("i8") * 1000 + T["j"][ix]) * 1000 + T["i"][ix]
    key = key * 8 + T["ott"][ix]
    u, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    s1 = np.bincount(inv, weights=e); s2 = np.bincount(inv, weights=e * e)
    ok = cnt >= 5
    mu = s1[ok] / cnt[ok]; var = (s2[ok] - cnt[ok] * mu ** 2) / (cnt[ok] - 1)
    nn = cnt[ok]
    sis = np.sum(nn * (mu ** 2 - var / nn)) / np.sum(nn)
    out["sistematico"] = float(np.sqrt(max(sis, 0)))
    out["sistematico_grezzo"] = float(np.sqrt(np.sum(nn * mu ** 2) / np.sum(nn)))
    for lo, hi, nome in ((0, 3, "sis0-3"), (3, 10, "sis3-10")):
        dm = np.bincount(inv, weights=d)[ok] / nn
        m = (dm >= lo) & (dm < hi)
        if m.sum() > 20:
            out[nome] = float(np.sqrt(max(np.sum(nn[m] * (mu[m] ** 2 - var[m] / nn[m])) / np.sum(nn[m]), 0)))
    # classi del rilievo (entro 10 km): scarto medio
    v = d < 10
    classi = {
        "ridosso alto (ang>4)": T["ang"][ix] > 4,
        "ridosso medio (1-4)": (T["ang"][ix] > 1) & (T["ang"][ix] <= 4),
        "da terra bassa": (T["terra"][ix] < 29.9) & (T["ang"][ix] <= 1),
        "dal mare": T["terra"][ix] >= 29.9,
        "valle/passo (>250 m)": T["valle"][ix] > 250,
        "canale (<15 km)": T["canale"][ix] < 15,
        "terra davanti <3 km": T["sotto"][ix] < 3,
        "sera (17 UTC)": T["sera"][ix] > 0,
        "mattina (05 UTC)": T["sera"][ix] == 0,
    }
    cb = {}
    for nome, m in classi.items():
        m = m & v
        if m.sum() > 500:
            cb[nome] = [float(np.mean(e[m])), int(m.sum())]
    out["classi"] = cb
    out["scarto_classi"] = float(np.mean([abs(x[0]) for x in cb.values()])) if cb else float("nan")
    return out


def main():
    tab, usc = sys.argv[1:3]
    args = sys.argv[3:]
    bers = "rel"
    if "--bersaglio" in args:
        bers = args[args.index("--bersaglio") + 1]
    cassa = "--cassaforte" in args
    nomi = [a for a in args if not a.startswith("--") and a not in ("rel", "ass")]
    os.makedirs(usc, exist_ok=True)
    T = carica(tab)
    fis = aggiungi_fisica(T)
    if fis:
        print("varianti fisiche:", fis)
    metodi = nomi[0].split(",") if nomi else list(METODI)
    y_all = T["y"] / (T["k0"] if bers == "rel" else 1.0)
    w_all = (T["mod"] * (T["k0"] if bers == "rel" else 1.0)) ** 2
    y_all = np.clip(y_all, 0.2, 3.0)
    usa = T["cassa"] if cassa else ~T["cassa"]
    risultati = {}
    oof = {}
    for nm in metodi:
        t0 = time.time()
        pred = np.full(len(T["y"]), np.nan, dtype="f4")
        per_area = {}
        for ia, area in enumerate(AREE):
            if cassa:
                tr = np.nonzero((T["area"] != ia) & ~T["cassa"])[0]    # allenato su tutto il resto del tempo
            else:
                tr = np.nonzero((T["area"] != ia) & usa)[0]
            te = np.nonzero((T["area"] == ia) & usa)[0]
            f = METODI[nm](T, tr, y_all[tr], w_all[tr])
            pred[te] = f(te)
            per_area[area] = misure(T, te, pred[te], bers)
        tutte = np.nonzero(usa)[0]
        risultati[nm] = {"tutte": misure(T, tutte, pred[tutte], bers), "aree": per_area, "secondi": time.time() - t0}
        oof[nm] = pred
        r = risultati[nm]["tutte"]
        print(f"{nm:22s} media {r['media']:+.2f} rmse {r['rmse']:.3f} sistematico {r['sistematico']:.3f} "
              f"(0-3 km {r.get('sis0-3', float('nan')):.3f}, 3-10 {r.get('sis3-10', float('nan')):.3f}) "
              f"classi {r['scarto_classi']:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    suff = ("_cassa" if cassa else "") + f"_{bers}"
    json.dump(risultati, open(os.path.join(usc, f"confronto{suff}.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(usc, f"previsioni{suff}.npz"), **oof)


if __name__ == "__main__":
    main()
