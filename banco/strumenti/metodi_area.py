"""Correzioni imparate dal radar area per area (si usano solo il modello e la scena: il radar serve per imparare).

  VERO=<radar_vero.npy> python3 metodi_area.py <tabella.npz> <uscita> [metodi] [--prova cassaforte|settimane] [k=v ...]

Metodi (tutti per area, per cella; g = fattore sul vento del modello letto dal mare):
  atlante     media pesata per cella e ottante della direzione del modello, tirata verso la media di area e ottante
  armoniche   per cella, g(θ) = c + a1 cos θ + b1 sin θ + a2 cos 2θ + b2 sin 2θ (θ = direzione del modello nella
              cella), minimi quadrati pesati con un freno sulle armoniche e un po' di lisciatura tra celle vicine
  analoghi    per ogni scena, le K scene passate più simili per vento al largo (forza e direzione), stabilità e ora;
              il fattore è la media del radar ÷ modello su quelle scene, tirata verso le armoniche dove sono poche
Prove:
  settimane   si impara sulle settimane pari fuori cassaforte e si giudica sulle dispari (per scegliere i parametri)
  cassaforte  si impara su tutto fuori cassaforte e si giudica sulle settimane della cassaforte (prova finale)
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confronto_metodi import AREE, KN, carica, misure

PAR = {"k_atl": 10.0, "lam1": 20.0, "lam2": 60.0, "liscia": 0.7, "K": 25, "k_an": 3.0, "sd_dir": 25.0,
       "sd_v": 2.5, "sd_st": 2.0, "sd_ora": 0.7, "vmin_an": 3.0, "sd_925": 0.0}


class Area:
    """Le scene e le celle di un'area in matrici scena × cella."""

    def __init__(self, T, ia, righe):
        self.T = T
        self.r = righe
        sc = T["scena"][righe]
        self.scene, inv = np.unique(sc, return_inverse=True)
        j = T["j"][righe].astype(int); i = T["i"][righe].astype(int)
        self.ny, self.nx = j.max() + 1, i.max() + 1
        self.cella = j * self.nx + i
        self.inv = inv
        ns, nc = len(self.scene), self.ny * self.nx
        self.R = np.full((ns, nc), np.nan, "f4"); self.M = np.full((ns, nc), np.nan, "f4"); self.D = np.full((ns, nc), np.nan, "f4")
        self.R[inv, self.cella] = T["radar"][righe]; self.M[inv, self.cella] = T["mod"][righe]; self.D[inv, self.cella] = T["dir"][righe]
        prima = np.unique(inv, return_index=True)[1]
        g = lambda k: T[k][righe][prima] if k in T else np.zeros(ns, "f4")
        self.t = g("t"); self.S0 = g("S0"); self.D0 = g("D0"); self.sera = g("sera")
        self.dTam = np.nan_to_num(g("dTam")); self.dT925 = np.nan_to_num(g("dT925"))
        self.cassa = g("cassa").astype(bool)
        self.v925 = np.nan_to_num(g("v925")); self.d925 = np.nan_to_num(g("d925"))
        self.sett = (self.t // (7 * 86400)).astype(int)

    def a_righe(self, G):
        """Da una matrice scena × cella di fattori alle righe della tabella."""
        return G[self.inv, self.cella]


def _media_ottanti(A, tr, p):
    R, M, D = A.R[tr], A.M[tr], A.D[tr]
    ok = np.isfinite(R) & np.isfinite(M)
    w = np.where(ok, M * M, 0.0); wy = np.where(ok, R * M, 0.0)
    ott = np.where(ok, np.round(np.nan_to_num(D) / 45.0).astype(int) % 8, 0)
    nc = R.shape[1]
    sw = np.zeros((8, nc)); swy = np.zeros((8, nc)); n = np.zeros((8, nc))
    for o in range(8):
        m = ok & (ott == o)
        sw[o] = np.where(m, w, 0).sum(0); swy[o] = np.where(m, wy, 0).sum(0); n[o] = m.sum(0)
    glob = swy.sum() / max(sw.sum(), 1e-9)
    so = sw.sum(1); syo = swy.sum(1); no = n.sum(1)
    wbo = so / np.maximum(no, 1)
    go = (syo + 30 * wbo * glob) / np.maximum(so + 30 * wbo, 1e-9)
    wbc = sw / np.maximum(n, 1)
    gc = (swy + p["k_atl"] * wbc * go[:, None]) / np.maximum(sw + p["k_atl"] * wbc, 1e-9)
    gc = np.where(n > 0, gc, go[:, None])
    return gc, glob


def atlante(A, tr, te, p):
    gc, _ = _media_ottanti(A, tr, p)
    D = A.D[te]
    ott = np.round(np.nan_to_num(D) / 45.0).astype(int) % 8
    return gc[ott, np.arange(gc.shape[1])[None, :]]


def _armoniche(A, tr, p):
    """Coefficienti per cella (5, nc): minimi quadrati pesati con freno, poi lisciatura spaziale leggera."""
    R, M, D = A.R[tr], A.M[tr], A.D[tr]
    ok = np.isfinite(R) & np.isfinite(M)
    th = np.radians(np.nan_to_num(D))
    w = np.where(ok, M * M, 0.0); y = np.where(ok, R / np.where(ok, M, 1), 0.0)
    y = np.clip(y, 0.2, 3.0)
    X = np.stack([np.ones_like(th), np.cos(th), np.sin(th), np.cos(2 * th), np.sin(2 * th)], axis=0)   # (5, ns, nc)
    glob = np.sum(w * y) / max(np.sum(w), 1e-9)
    wbar = np.sum(w) / max(ok.sum(), 1)
    pen = np.array([1e-3, p["lam1"], p["lam1"], p["lam2"], p["lam2"]]) * wbar
    nc = R.shape[1]
    XtWX = np.einsum("asc,bsc,sc->cab", X, X, w)
    XtWy = np.einsum("asc,sc->ca", X, w * (y - glob))
    XtWX[:, range(5), range(5)] += pen[None, :]
    # il termine costante è tirato verso glob da un freno pari a k_atl scene medie
    XtWX[:, 0, 0] += p["k_atl"] * wbar
    coef = np.linalg.solve(XtWX, XtWy[..., None])[..., 0].T          # (5, nc)
    coef[0] += glob
    if p["liscia"] > 0:
        C = coef.reshape(5, A.ny, A.nx)
        coef = np.stack([gaussian_filter(C[k], p["liscia"], mode="nearest") for k in range(5)]).reshape(5, nc)
    return coef


def armoniche(A, tr, te, p, coef=None):
    if coef is None:
        coef = _armoniche(A, tr, p)
    th = np.radians(np.nan_to_num(A.D[te]))
    X = np.stack([np.ones_like(th), np.cos(th), np.sin(th), np.cos(2 * th), np.sin(2 * th)], axis=0)
    return np.einsum("asc,ac->sc", X, coef)


def _vett_scena(A, ix, p):
    a = np.radians(A.D0[ix]); b = np.radians(A.d925[ix])
    col = [A.S0[ix] * np.sin(a) / p["sd_v"], A.S0[ix] * np.cos(a) / p["sd_v"],
           A.dT925[ix] / p["sd_st"], A.dTam[ix] / p["sd_st"], A.sera[ix] / p["sd_ora"]]
    if p["sd_925"] > 0:
        col += [A.v925[ix] * np.sin(b) / p["sd_925"], A.v925[ix] * np.cos(b) / p["sd_925"]]
    return np.stack(col, axis=1)


def analoghi(A, tr, te, p):
    coef = _armoniche(A, tr, p)
    base = armoniche(A, tr, te, p, coef)
    Xtr = _vett_scena(A, tr, p); Xte = _vett_scena(A, te, p)
    R, M = A.R[tr], A.M[tr]
    ok = np.isfinite(R) & np.isfinite(M) & (M > p["vmin_an"])
    wy = np.where(ok, R * M, 0.0); ww = np.where(ok, M * M, 0.0)
    G = np.empty_like(base)
    K = int(p["K"])
    wbar = np.sum(ww) / max(ok.sum(), 1)
    for s in range(len(te)):
        d2 = np.sum((Xtr - Xte[s]) ** 2, axis=1)
        idx = np.argpartition(d2, min(K, len(d2) - 1))[:K]
        pk = np.exp(-0.5 * d2[idx] / max(np.median(d2[idx]), 1e-6))
        num = (pk[:, None] * wy[idx]).sum(0); den = (pk[:, None] * ww[idx]).sum(0)
        lam = p["k_an"] * wbar
        G[s] = (num + lam * base[s]) / (den + lam)
    return G


METODI = {"modello": lambda A, tr, te, p: np.ones((len(te), A.R.shape[1])),
          "atlante": atlante, "armoniche": armoniche, "analoghi": analoghi}


def main():
    tab, usc = sys.argv[1:3]
    args = sys.argv[3:]
    prova = args[args.index("--prova") + 1] if "--prova" in args else "settimane"
    p = dict(PAR)
    for a in args:
        if "=" in a:
            k, v = a.split("=", 1); p[k] = float(v)
    nomi = [a for a in args if "=" not in a and not a.startswith("--") and a not in ("settimane", "cassaforte")]
    metodi = nomi[0].split(",") if nomi else list(METODI)
    os.makedirs(usc, exist_ok=True)
    T = carica(tab)
    aree_in = sorted(set(T["area"].tolist()))
    pred = {m: np.full(len(T["radar"]), np.nan, "f4") for m in metodi}
    giudicate = np.zeros(len(T["radar"]), bool)
    t0 = time.time()
    for ia in aree_in:
        righe = np.nonzero(T["area"] == ia)[0]
        A = Area(T, ia, righe)
        if prova == "cassaforte":
            tr = np.nonzero(~A.cassa)[0]; te = np.nonzero(A.cassa)[0]
        else:
            tr = np.nonzero(~A.cassa & (A.sett % 2 == 0))[0]; te = np.nonzero(~A.cassa & (A.sett % 2 == 1))[0]
        sel_te = np.isin(A.inv, te)
        pos = np.full(len(A.scene), -1); pos[te] = np.arange(len(te))
        for m in metodi:
            G = METODI[m](A, tr, te, p)
            r = righe[sel_te]
            pred[m][r] = G[pos[A.inv[sel_te]], A.cella[sel_te]]
        giudicate[righe[sel_te]] = True
    ix = np.nonzero(giudicate)[0]
    ris = {"prova": prova, "parametri": p}
    for m in metodi:
        r = misure(T, ix, pred[m][ix], "ass")
        per_area = {}
        for ia in aree_in:
            jx = ix[T["area"][ix] == ia]
            if len(jx):
                per_area[AREE[ia]] = misure(T, jx, pred[m][jx], "ass")
        ris[m] = {"tutte": r, "aree": per_area}
        print(f"{m:12s} rmse {r['rmse']:.3f} sistematico {r['sistematico']:.3f} (0-3 km {r.get('sis0-3', float('nan')):.3f}, "
              f"3-10 {r.get('sis3-10', float('nan')):.3f}) media {r['media']:+.2f}", flush=True)
    nome = f"area_{prova}_" + "_".join(f"{k}{v:g}" for k, v in sorted(p.items()) if v != PAR.get(k))
    json.dump(ris, open(os.path.join(usc, nome.rstrip("_") + ".json"), "w"), indent=1)
    print(f"[{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
