"""Per area e direzione di provenienza del modello (letto dal mare): radar ÷ modello vicino alla costa.
Per non cadere nel tranello della regressione verso la media (condizionare sulla forza del modello stesso
gonfia il rapporto a vento debole e lo sgonfia a vento forte), la forza si prende dal largo: la mediana del modello
nelle celle oltre 10 km della stessa scena.

  python3 diagnosi_direzioni.py <cartella del metro> area [area ...]
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diagnosi_rilievo import campi
KN = 1.943844
SETT = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO"]
cart = sys.argv[1]
for area in sys.argv[2:]:
    z, mu, mv = campi(cart, area)
    R = z["radar"].astype("f8"); M = np.hypot(mu, mv); D = np.broadcast_to(z["dist"], R.shape)
    larg = np.where(np.isfinite(M) & (D > 10), M, np.nan)
    Ulargo = np.nanmedian(larg.reshape(len(R), -1), axis=1)            # forza al largo, per scena
    ulm = np.where(np.isfinite(M) & (D > 10), mu, np.nan).reshape(len(R), -1); vlm = np.where(np.isfinite(M) & (D > 10), mv, np.nan).reshape(len(R), -1)
    dirs = (np.degrees(np.arctan2(-np.nanmedian(ulm, axis=1), -np.nanmedian(vlm, axis=1))) + 360) % 360
    k = np.round(dirs / 22.5).astype(int) % 16
    print(f"\n== {area}: per direzione del modello al largo (scene con oltre 8 nodi al largo)")
    print(f"   {'da':<4} {'scene':>5}  {'0-3 km':>16}  {'3-10 km':>16}  {'oltre 10':>16}   (radar/modello, errore nodi)")
    for s in range(16):
        sc = (k == s) & (Ulargo * KN > 8)
        if sc.sum() < 3:
            continue
        cols = []
        for lo, hi in [(0, 3), (3, 10), (10, 99)]:
            m = sc[:, None, None] & np.isfinite(R) & np.isfinite(M) & (M > 1.5) & (D >= lo) & (D < hi)
            cols.append(f"{np.median(R[m] / M[m]):5.2f} {np.mean(M[m] - R[m]) * KN:+5.1f}" if m.sum() > 50 else "      -     ")
        print(f"   {SETT[s]:<4} {sc.sum():5d}  " + "  ".join(f"{c:>16}" for c in cols))
