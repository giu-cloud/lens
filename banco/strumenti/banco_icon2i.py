"""ICON-2I sulle celle del radar: lo stesso confronto fatto con ECMWF, con il modello fine (storia di Open-Meteo).

  python3 banco_icon2i.py <tabella.npz> <cartella icon2i_*.npz> <cartella del metro> <rilievo .rgr> <uscita.npz> [--tutti]

Per ogni passaggio del radar con ICON-2I disponibile (dall'aprile 2025), il vento di ICON-2I si porta sulle celle:
  - nel tempo, lineare tra le due ore piene vicine;
  - nello spazio, bilineare dai soli punti di mare della griglia di ICON-2I (0,02° × 0,025°): i punti con più di metà
    terra (rilievo di Shift) non si usano, come fa Shift con gli altri modelli; se i quattro sono a terra, il punto di
    mare più vicino entro due maglie.
Scrive una tabella con le stesse righe della tabella di ECMWF che hanno ICON-2I: le colonne del modello (mod, dir, S0, D0,
dDir, rS0) prendono i valori di ICON-2I; quelle di ECMWF restano come mod_ecmwf, dir_ecmwf.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confronto_metodi import AREE
from rilievo import Rilievo


def terra_nodi(ril, lat, lon):
    s, n, w, e = lat.min() - 0.05, lat.max() + 0.05, lon.min() - 0.05, lon.max() + 0.05
    rl, ro, H = ril.griglia(0, s, n, w, e)
    T = np.zeros((len(lat), len(lon)), bool)
    c = rl[1] - rl[0]
    for a, la in enumerate(lat):
        ja = int((la - 0.01 - rl[0]) / c); jb = int((la + 0.01 - rl[0]) / c) + 1
        for b, lo in enumerate(lon):
            ia = int((lo - 0.0125 - ro[0]) / c); ib = int((lo + 0.0125 - ro[0]) / c) + 1
            T[a, b] = (H[max(ja, 0):jb, max(ia, 0):ib] > 0).mean() > 0.5
    return T


def pesi_mare(lat, lon, terra, lat_t, lon_t):
    """Per ogni cella di arrivo: indici e pesi (fino a 4) dei nodi di mare."""
    dla = lat[1] - lat[0]; dlo = lon[1] - lon[0]
    idx = np.zeros((len(lat_t), len(lon_t), 4, 2), int); w = np.zeros((len(lat_t), len(lon_t), 4))
    for a, la in enumerate(lat_t):
        fj = (la - lat[0]) / dla; j0 = int(np.floor(fj)); y = fj - j0
        for b, lo in enumerate(lon_t):
            fi = (lo - lon[0]) / dlo; i0 = int(np.floor(fi)); x = fi - i0
            pts = [(j0, i0, (1 - x) * (1 - y)), (j0, i0 + 1, x * (1 - y)), (j0 + 1, i0, (1 - x) * y), (j0 + 1, i0 + 1, x * y)]
            k = 0
            for j, i, ww in pts:
                if 0 <= j < len(lat) and 0 <= i < len(lon) and not terra[j, i] and ww > 0:
                    idx[a, b, k] = (j, i); w[a, b, k] = ww; k += 1
            if w[a, b].sum() > 1e-6:
                w[a, b] /= w[a, b].sum()
                continue
            best = None; bd = 1e9
            for j in range(j0 - 2, j0 + 4):
                for i in range(i0 - 2, i0 + 4):
                    if 0 <= j < len(lat) and 0 <= i < len(lon) and not terra[j, i]:
                        dd = (j - fj) ** 2 + (i - fi) ** 2
                        if dd < bd:
                            bd, best = dd, (j, i)
            if best:
                idx[a, b, 0] = best; w[a, b, 0] = 1.0
    return idx, w


def main():
    tab, cart_i, metro, ril_p, usc = sys.argv[1:6]
    T = dict(np.load(tab))
    ril = Rilievo(ril_p)
    n = len(T["radar"])
    I_mod = np.full(n, np.nan, "f4"); I_dir = np.full(n, np.nan, "f4"); I_S0 = np.full(n, np.nan, "f4"); I_D0 = np.full(n, np.nan, "f4")
    for ia, area in enumerate(AREE):
        p = os.path.join(cart_i, f"icon2i_{area}.npz")
        if not os.path.exists(p) or not (T["area"] == ia).any():
            continue
        Z = np.load(p)
        lat, lon, ti = Z["lat"], Z["lon"], Z["t"]
        o = np.argsort(ti); ti = ti[o]
        U = Z["u"][o].astype("f4"); V = Z["v"][o].astype("f4")
        z = np.load(os.path.join(metro, f"campi_{area}.npz"))
        terra = terra_nodi(ril, lat, lon)
        idx, w = pesi_mare(lat, lon, terra, z["lat"], z["lon"])
        dist = z["dist"]
        righe = np.nonzero(T["area"] == ia)[0]
        sc = T["scena"][righe]
        scene = np.unique(sc)
        fatte = 0
        for s in scene:
            rr = righe[sc == s]
            t = T["t"][rr[0]]
            k = np.searchsorted(ti, t)
            if k == 0 or k >= len(ti) or ti[k] - ti[k - 1] > 3600 or t - ti[k - 1] > 3600:
                continue
            f = (t - ti[k - 1]) / (ti[k] - ti[k - 1])
            u = (1 - f) * U[k - 1] + f * U[k]; v = (1 - f) * V[k - 1] + f * V[k]
            uc = np.sum(u[idx[..., 0], idx[..., 1]] * w, axis=-1); vc = np.sum(v[idx[..., 0], idx[..., 1]] * w, axis=-1)
            j = T["j"][rr].astype(int); i = T["i"][rr].astype(int)
            I_mod[rr] = np.hypot(uc[j, i], vc[j, i])
            I_dir[rr] = (np.degrees(np.arctan2(-uc[j, i], -vc[j, i])) + 360) % 360
            largo = dist >= 8 if (dist >= 8).sum() >= 20 else dist >= np.quantile(dist[dist > 0], 0.75)
            ul, vl = uc[largo].mean(), vc[largo].mean()
            I_S0[rr] = np.median(np.hypot(uc, vc)[largo]); I_D0[rr] = (np.degrees(np.arctan2(-ul, -vl)) + 360) % 360
            fatte += 1
        print(f"{area}: {fatte} passaggi con ICON-2I su {len(scene)}", flush=True)
    ok = np.isfinite(I_mod) & (I_mod > 2.0)
    out = {k: v[ok] for k, v in T.items()}
    out["mod_ecmwf"] = T["mod"][ok]; out["dir_ecmwf"] = T["dir"][ok]; out["S0_ecmwf"] = T["S0"][ok]; out["D0_ecmwf"] = T["D0"][ok]
    out["mod"] = I_mod[ok]; out["dir"] = I_dir[ok]; out["S0"] = I_S0[ok]; out["D0"] = I_D0[ok]
    out["dDir"] = ((out["dir"] - out["D0"] + 540) % 360 - 180).astype("f4")
    out["riga"] = np.nonzero(ok)[0]
    np.savez(usc, **out)
    print(f"{ok.sum()} righe con ICON-2I -> {usc}")


if __name__ == "__main__":
    main()
