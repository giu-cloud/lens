"""Prova incrociata di una correzione locale (fattore per cella e settore di direzione) sopra ECMWF letto dal mare.
Metà dei passaggi (giorni pari) fa le mappe, l'altra metà le prova, e viceversa."""
import sys
import numpy as np
from scipy import ndimage

KN = 1.943844
BINS = [(0, 3), (3, 10), (10, 999)]


def mappe(R, M, U, V, idx, n_sett, sig, k):
    """Fattore per settore: somma radar / somma modello, lisciato, tirato verso 1 con k passaggi a priori."""
    out = {}
    d = (np.degrees(np.arctan2(-U[idx], -V[idx])) + 360) % 360
    s = ((d + 180 / n_sett) // (360 / n_sett)).astype(int) % n_sett
    for q in range(n_sett):
        ii = idx[s == q]
        if len(ii) == 0:
            continue
        ok = np.isfinite(R[ii]) & np.isfinite(M[ii]) & (M[ii] > 1.5)
        num = np.where(ok, R[ii], 0).sum(0); den = np.where(ok, M[ii], 0).sum(0); cnt = ok.sum(0)
        f = np.where(den > 0, num / np.maximum(den, 1e-6), 1.0)
        w = cnt / (cnt + k)
        f = 1 + (f - 1) * w
        if sig > 0:
            f = ndimage.gaussian_filter(np.nan_to_num(f, nan=1.0), sig)
        out[q] = np.clip(f, 0.6, 1.8)
    return out


def main():
    base = sys.argv[1]
    for area in sys.argv[2:]:
        z = np.load(f"{base}/campi_{area}.npz")
        R = z["radar"].astype("f8"); M = z["mare"].astype("f8"); dist = z["dist"]; t = z["t"]
        # direzione del modello per scena: dal file di campi non abbiamo u, v; si rifà dal settore del radar medio? no: dal cmems
        # (u, v del modello non salvati): uso la direzione del vento ECMWF del prodotto dalla scena originale
        import costa_analisi as ca
        lat, lon, scene, _ = ca.carica("dati", area)
        tt = {s["t"]: s for s in scene}
        U = np.array([np.nanmean(tt[x]["mu"]) for x in t]); V = np.array([np.nanmean(tt[x]["mv"]) for x in t])
        giorno = (t // 86400) % 2
        D = np.broadcast_to(dist, R.shape)
        for n_sett, sig, k in [(8, 1.5, 5), (8, 3, 5), (4, 3, 5), (8, 3, 15), (16, 3, 10)]:
            res = np.zeros((len(BINS), 2, 2))      # bin, (prima, dopo), (somma e^2, n)
            for test in (0, 1):
                tr = np.where(giorno != test)[0]; te = np.where(giorno == test)[0]
                F = mappe(R, M, U, V, tr, n_sett, sig, k)
                d = (np.degrees(np.arctan2(-U[te], -V[te])) + 360) % 360
                s = ((d + 180 / n_sett) // (360 / n_sett)).astype(int) % n_sett
                for j, kk in enumerate(te):
                    f = F.get(s[j], np.ones_like(dist))
                    ok = np.isfinite(R[kk]) & np.isfinite(M[kk]) & (M[kk] > 1.5)
                    for b, (lo, hi) in enumerate(BINS):
                        m = ok & (dist >= lo) & (dist < hi)
                        res[b, 0, 0] += ((M[kk][m] - R[kk][m]) ** 2).sum(); res[b, 0, 1] += m.sum()
                        res[b, 1, 0] += ((M[kk][m] * f[m] - R[kk][m]) ** 2).sum(); res[b, 1, 1] += m.sum()
            print(f"{area:<20} settori {n_sett:2d} liscio {sig} k {k:2d}: " + "  ".join(
                f"{lo}-{hi} km {np.sqrt(res[b, 0, 0] / res[b, 0, 1]) * KN:.2f} -> {np.sqrt(res[b, 1, 0] / res[b, 1, 1]) * KN:.2f}" for b, (lo, hi) in enumerate(BINS)))


if __name__ == "__main__":
    main()
