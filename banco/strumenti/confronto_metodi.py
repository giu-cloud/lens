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
NSUB = int(os.environ.get("NSUB", "500000"))


def carica(p):
    T = dict(np.load(p))
    vero = os.environ.get("VERO")
    if vero:
        # bersaglio: il vento vero a 10 m (il radar convertito con la stabilità di ECMWF, vento_vero.py), non quello neutro
        T["radar_neutro"] = T["radar"]
        T["radar"] = np.load(vero)
    T["y"] = (T["radar"] / T["mod"]).astype("f4")
    for k in [k for k in T if k.startswith("q_") or k.startswith("sw_")]:
        T[k] = np.nan_to_num(T[k], nan=1.0)
    T["lmin"] = np.minimum(T["lsx"], T["ldx"])
    T["canale"] = np.minimum(T["lsx"] + T["ldx"], 60.0)
    T["sera"] = (T["ora"] >= 12).astype("f4")
    T["msin"] = np.sin(2 * np.pi * T["mese"] / 12).astype("f4")
    T["mcos"] = np.cos(2 * np.pi * T["mese"] / 12).astype("f4")
    T["rS0"] = (T["mod"] / np.maximum(T["S0"], 0.5)).astype("f4")
    T["ott"] = (np.round(T["dir"] / 45.0).astype(int) % 8).astype("i1")
    T["cassa"] = ((T["t"] // SETTIMANA) % 5 == 2)
    if "d925" in T:
        T["d925r"] = ((T["d925"] - T["D0"] + 540) % 360 - 180).astype("f4")     # rotazione tra 925 hPa e superficie
    return T


STATICHE = ["dist", "terra", "h15", "ang", "hc", "valle", "f5", "f15", "sotto", "hsotto", "angsotto", "lmin", "canale",
            "angw", "fw15", "h10", "hm10", "terra10"]
SCENA = ["S0", "rS0", "dDir", "sera", "msin", "mcos"]
STAB = ["sst", "dTam", "dT925", "dT850", "terra_mare", "taglio", "raff", "v925", "d925r", "v850", "umid"]


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
        sub = tr if len(tr) <= NSUB else np.sort(rng.choice(tr, NSUB, replace=False))
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


def _potenza(col):
    """g = c × q^a, con q una risposta fisica (Lens, strato marino): due numeri soli, allenati sulle altre aree."""
    def fit(T, tr, y, w):
        q = np.clip(np.nan_to_num(T[col][tr], nan=1.0), 0.2, 5.0)
        lq = np.log(q)
        best = None
        for a in np.linspace(-0.5, 1.5, 41):
            p = np.exp(a * lq)
            c = np.sum(w * y * p) / np.sum(w * p * p)
            r = np.sum(w * (y - c * p) ** 2)
            if best is None or r < best[0]:
                best = (r, a, c)
        _, a, c = best
        print(f"    {col}: a={a:.2f} c={c:.3f}", flush=True) if os.environ.get("VERBOSO") else None
        return lambda te: (c * np.clip(np.nan_to_num(T[col][te], nan=1.0), 0.2, 5.0) ** a).astype("f4")
    return fit


def _combinata(cols, extra=()):
    """log g = c + Σ a_i log q_i (+ termini extra lineari): i moduli fisici pesati, pochi numeri, Gauss-Newton sul
    rapporto pesato (lo stesso obiettivo degli altri metodi)."""
    def X_(T, ix):
        X = [np.ones(len(ix))]
        for c in cols:
            X.append(np.log(np.clip(np.nan_to_num(T[c][ix], nan=1.0), 0.2, 5.0)))
        for c in extra:
            X.append(np.nan_to_num(T[c][ix], nan=0.0))
        return np.stack(X, axis=1)

    def fit(T, tr, y, w):
        X = X_(T, tr)
        b = np.zeros(X.shape[1]); b[0] = np.log(np.sum(w * y) / np.sum(w))
        lam = 1e-4 * np.sum(w)
        for _ in range(12):
            p = np.exp(X @ b)
            J = X * p[:, None]
            A = (J * w[:, None]).T @ J + lam * np.eye(len(b)); A[0, 0] -= lam
            db = np.linalg.solve(A, (J * w[:, None]).T @ (y - p))
            b = b + db
            if np.max(np.abs(db)) < 1e-5:
                break
        if os.environ.get("VERBOSO"):
            print("    combinata:", dict(zip(["c"] + list(cols) + list(extra), np.round(b, 3))), flush=True)
        return lambda te: np.exp(X_(T, te) @ b).astype("f4")
    return fit


def m_atlante(T, tr, y, w, k_shr=10.0, blocco=1):
    """Carta imparata dal radar per area, cella e ottante (solo prova nel tempo): rapporto pesato, tirato verso la media
    di area e ottante quando i passaggi sono pochi (k_shr passaggi equivalenti), e quella verso la media generale."""
    def chiavi(ix):
        a = T["area"][ix].astype("i8"); o = T["ott"][ix].astype("i8")
        jj = T["j"][ix].astype("i8") // blocco; ii = T["i"][ix].astype("i8") // blocco
        return ((a * 1000 + jj) * 1000 + ii) * 8 + o, a * 8 + o
    kc, ko = chiavi(tr)
    glob = np.sum(w * y) / np.sum(w)
    def media(k):
        u, inv = np.unique(k, return_inverse=True)
        sw = np.bincount(inv, weights=w); swy = np.bincount(inv, weights=w * y); n = np.bincount(inv)
        return u, sw, swy, n
    uo, swo, swyo, no = media(ko)
    wbo = swo / np.maximum(no, 1)
    go = (swyo + 30 * wbo * glob) / (swo + 30 * wbo)
    uc, swc, swyc, nc = media(kc)
    # la media di area e ottante per ogni cella
    pos = np.searchsorted(uo, (uc // (8 * 10 ** 6)) * 8 + uc % 8)      # la chiave di area e ottante della cella
    pos = np.clip(pos, 0, len(uo) - 1)
    base = go[pos]
    wbc = swc / np.maximum(nc, 1)
    gc = (swyc + k_shr * wbc * base) / (swc + k_shr * wbc)
    def predici(te):
        k1, k2 = chiavi(te)
        p1 = np.searchsorted(uc, k1); p1 = np.clip(p1, 0, len(uc) - 1)
        p2 = np.searchsorted(uo, k2); p2 = np.clip(p2, 0, len(uo) - 1)
        g = np.where(uc[p1] == k1, gc[p1], np.where(uo[p2] == k2, go[p2], glob))
        return g.astype("f4")
    return predici


METODI = {}


def costruisci_metodi(T):
    """I metodi, con le caratteristiche che la tabella ha."""
    st = STATICHE + [k for k in ("cresta", "passo", "tsx", "tdx") if k in T]
    sc = SCENA + [k for k in STAB if k in T and np.isfinite(T[k]).mean() > 0.9]
    fis = [k for k in T if k.startswith("q_") and k != "q_lens054" and np.isfinite(T[k]).all()]
    sw = [k for k in ("sw_loc", "sw_lett", "sw_spess") if k in T and np.isfinite(T[k]).all()]
    METODI.update({
        "atlante": m_atlante,
        "atlante_k30": lambda T_, tr, y, w: m_atlante(T_, tr, y, w, k_shr=30.0),
        "atlante_b2": lambda T_, tr, y, w: m_atlante(T_, tr, y, w, k_shr=10.0, blocco=2),
        "atlante_b3k20": lambda T_, tr, y, w: m_atlante(T_, tr, y, w, k_shr=20.0, blocco=3),
        "modello": m_modello, "lens054": m_lens054, "costante": m_costante, "distanza": m_distanza,
        "dist_rilievo": m_dist_rilievo, "lineare": m_lineare,
        "alberi_rilievo": _hgb(st),
        "alberi_tutto": _hgb(st + SCENA),
        "alberi_tutto_lens": _hgb(st + SCENA + ["q_lens054"]),
    })
    for k in fis:
        METODI["potenza_" + k[2:]] = _potenza(k)
        METODI["diretto_" + k[2:]] = (lambda col: (lambda T_, tr, y, w: (lambda te: T_[col][te].astype("f4"))))(k)
    fis_ok = [k for k in ["q_lens054"] + fis if k in T]
    METODI["combinata"] = _combinata(fis_ok)
    if fis:
        METODI["alberi_rilievo_fisica"] = _hgb(st + fis + sw)
        METODI["alberi_tutto_fisica"] = _hgb(st + SCENA + fis + sw)
    if len(sc) > len(SCENA):
        METODI["alberi_tutto_stab"] = _hgb(st + sc)
        METODI["alberi_tutto_stab_fisica"] = _hgb(st + sc + fis + sw)
    return fis + sw, sc


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
    tempo = "--tempo" in args
    nomi = [a for a in args if not a.startswith("--") and a not in ("rel", "ass")]
    os.makedirs(usc, exist_ok=True)
    T = carica(tab)
    fis, sc = costruisci_metodi(T)
    print("fisica:", fis, "scena:", sc)
    metodi = nomi[0].split(",") if nomi else list(METODI)
    y_all = T["y"] / (T["k0"] if bers == "rel" else 1.0)
    w_all = (T["mod"] * (T["k0"] if bers == "rel" else 1.0)) ** 2
    y_all = np.clip(y_all, 0.2, 3.0)
    usa = T["cassa"] if (cassa or tempo) else ~T["cassa"]
    risultati = {}
    oof = {}
    for nm in metodi:
        t0 = time.time()
        pred = np.full(len(T["y"]), np.nan, dtype="f4")
        per_area = {}
        f_tempo = None
        if tempo:
            tr = np.nonzero(~T["cassa"])[0]
            f_tempo = METODI[nm](T, tr, y_all[tr], w_all[tr])
        for ia, area in enumerate(AREE):
            if tempo:
                # prova nel tempo: allenato su tutte le aree fuori dalla cassaforte, giudicato sulle settimane della cassaforte
                tr = np.nonzero(~T["cassa"])[0]
            elif cassa:
                tr = np.nonzero((T["area"] != ia) & ~T["cassa"])[0]    # allenato su tutto il resto del tempo
            else:
                tr = np.nonzero((T["area"] != ia) & usa)[0]
            te = np.nonzero((T["area"] == ia) & usa)[0]
            f = f_tempo if f_tempo is not None else METODI[nm](T, tr, y_all[tr], w_all[tr])
            pred[te] = f(te)
            per_area[area] = misure(T, te, pred[te], bers)
        tutte = np.nonzero(usa)[0]
        risultati[nm] = {"tutte": misure(T, tutte, pred[tutte], bers), "aree": per_area, "secondi": time.time() - t0}
        oof[nm] = pred
        r = risultati[nm]["tutte"]
        print(f"{nm:22s} media {r['media']:+.2f} rmse {r['rmse']:.3f} sistematico {r['sistematico']:.3f} "
              f"(0-3 km {r.get('sis0-3', float('nan')):.3f}, 3-10 {r.get('sis3-10', float('nan')):.3f}) "
              f"classi {r['scarto_classi']:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    suff = ("_tempo" if tempo else "_cassa" if cassa else "") + f"_{bers}"
    json.dump(risultati, open(os.path.join(usc, f"confronto{suff}.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(usc, f"previsioni{suff}.npz"), **oof)


if __name__ == "__main__":
    main()
