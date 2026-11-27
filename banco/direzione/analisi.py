"""Analisi della direzione: errore assoluto medio (gradi) contro le stazioni, per versione di Lens, punto e condizione."""
import glob, sys, numpy as np
L = "/home/claude/lavoro/dir"
NOMI = ["M", "L0", "LA", "LM", "LS"]


def carica(suf=""):
    R = []
    for f in sorted(glob.glob(f"{L}/dir_*{suf}.npz")):
        a = f.split("dir_")[-1][:-4]
        if (suf == "") != ("_" not in a.replace("-", "")):
            continue
        Z = np.load(f); X = Z["X"]; st = Z["st"]
        for k in range(len(st)):
            R.append((a, st[k], X[k]))
    return R


def dire(u, v):
    return np.degrees(np.arctan2(-u, -v)) % 360


def diff(a, b):
    return (a - b + 180) % 360 - 180


def tabella(R, vmin=5.0, dove="c", solo_cassa=None, per="st"):
    A = np.array([r[2] for r in R]); st = np.array([r[1] for r in R]); ar = np.array([r[0] for r in R])
    t = A[:, 0]; vs = A[:, 1]; ds = A[:, 2]
    o = 5 if dove == "c" else 15
    UV = {n: (A[:, o + 2 * k], A[:, o + 2 * k + 1]) for k, n in enumerate(NOMI)}
    sm = np.hypot(*UV["M"])
    ok = (sm >= vmin) & (vs >= 2.0) & np.isfinite(UV["LA"][0])
    if solo_cassa is not None:
        ok &= ((t // (7 * 86400)) % 5 == 2) == solo_cassa
    chiavi = st if per == "st" else ar
    out = []
    for k in list(np.unique(chiavi[ok])) + ["TUTTE"]:
        m = ok & ((chiavi == k) if k != "TUTTE" else True)
        if m.sum() < 30:
            continue
        e = {n: diff(dire(*UV[n])[m], ds[m]) for n in NOMI}
        out.append((k, ar[m][0] if k != "TUTTE" else "", int(m.sum()), {n: (float(np.mean(np.abs(e[n]))), float(np.mean(e[n]))) for n in NOMI}))
    return out


if __name__ == "__main__":
    R = carica(sys.argv[1] if len(sys.argv) > 1 else "")
    for dove in ("c", "m"):
        for vmin in (5.0, 2.5):
            print(f"\n== {'cella della stazione' if dove == 'c' else 'mare aperto vicino'}, modello >= {vmin} m/s: errore assoluto medio in gradi (medio con segno)")
            for k, a, n, e in tabella(R, vmin, dove):
                print(f"{a:20s} {k:16s} {n:5d}  " + "  ".join(f"{x}: {e[x][0]:5.1f} ({e[x][1]:+5.1f})" for x in NOMI))
