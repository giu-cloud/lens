"""Diagnosi: dove e perché il radar vede più vento del modello letto dal mare, secondo il rilievo sopravvento.

  python3 diagnosi_rilievo.py <cartella del metro (campi_*.npz, sho_*.bin)> <rilievo di Shift .rgr> <uscita>

Per ogni cella del radar e per 16 direzioni di provenienza si cammina sopravvento fino a 30 km, sul rilievo di Shift
(464 m), e si misurano:
  terra   distanza dalla prima terra sopravvento (km; 30 = mare aperto)
  h15     quota massima sopravvento entro 15 km dalla costa (m)
  ang     angolo massimo sotto cui si vede il rilievo sopravvento dalla cella (gradi, come Winstral)
  hc      quota media dei primi 3 km di terra sopravvento (m)
  valle   quanto il rilievo ai lati (±4 km) sta sopra quello sulla linea del vento, 2-6 km dentro terra (m):
          positivo = il vento arriva da una valle o da un passo
Poi per ogni scena si prende la direzione del modello letto dal mare e si confrontano radar e modello.
Il rapporto è radar ÷ modello (mediana), l'errore è modello − radar in nodi.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rilievo import Rilievo

KN = 1.943844
AREE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]
NSET = 16
PASSO = 0.25                      # km
RMAX = 30.0


def campi(cart, area):
    z = np.load(os.path.join(cart, f"campi_{area}.npz"))
    raw = np.fromfile(os.path.join(cart, f"sho_{area}.bin"), dtype="<f4", offset=16)
    ny, nx = len(z["lat"]), len(z["lon"])
    raw = raw.reshape(-1, 8, ny, nx)
    return z, raw[:, 2].astype("f8"), raw[:, 3].astype("f8")


def caratteristiche(lat, lon, ril):
    """Per ogni cella (ny, nx) e direzione di provenienza (16): terra, h15, ang, hc, valle."""
    s, n, w, e = lat.min() - 0.35, lat.max() + 0.35, lon.min() - 0.45, lon.max() + 0.45
    rl, ro, H = ril.griglia(0, s, n, w, e)
    H = np.maximum(H, 0).astype("f8")
    c = rl[1] - rl[0]
    kmy = 111.19; kmx = 111.19 * np.cos(np.radians(lat.mean()))
    LA, LO = np.meshgrid(lat, lon, indexing="ij")

    def quota(la, lo):
        j = np.clip(np.round((la - rl[0]) / c).astype(int), 0, len(rl) - 1)
        i = np.clip(np.round((lo - ro[0]) / c).astype(int), 0, len(ro) - 1)
        return H[j, i]

    passi = np.arange(PASSO, RMAX + 1e-9, PASSO)
    out = {k: np.zeros((NSET,) + LA.shape) for k in ("terra", "h15", "ang", "hc", "valle")}
    for k in range(NSET):
        a = np.radians(k * 360.0 / NSET)                  # direzione da cui viene il vento
        dy, dx = np.cos(a), np.sin(a)                     # verso sopravvento
        py, px = -dx, dy                                  # perpendicolare
        hr = np.stack([quota(LA + r * dy / kmy, LO + r * dx / kmx) for r in passi])      # (nr, ny, nx)
        terra = hr > 0
        prima = np.where(terra.any(axis=0), terra.argmax(axis=0), len(passi))
        dterra = np.where(prima < len(passi), passi[np.minimum(prima, len(passi) - 1)], RMAX)
        r3 = passi[:, None, None]
        dentro15 = (r3 >= dterra[None]) & (r3 <= dterra[None] + 15.0)
        h15 = np.where(dentro15, hr, 0).max(axis=0)
        ang = np.degrees(np.arctan(hr / (r3 * 1000.0))).max(axis=0)
        dentro3 = (r3 >= dterra[None]) & (r3 <= dterra[None] + 3.0)
        hc = np.where(dentro3, hr, 0).sum(axis=0) / np.maximum(dentro3.sum(axis=0), 1)
        # valle: 2-6 km dentro terra, rilievo sulla linea contro il massimo ai lati (±4 km)
        val = np.zeros(LA.shape); nval = np.zeros(LA.shape)
        for dd in (2.0, 4.0, 6.0):
            r = np.minimum(dterra + dd, RMAX)
            cy = LA + r * dy / kmy; cx = LO + r * dx / kmx
            linea = quota(cy, cx)
            lati = np.max([quota(cy + q * py / kmy, cx + q * px / kmx) for q in (-4, -3, -2, 2, 3, 4)], axis=0)
            ok = dterra < RMAX
            val += np.where(ok, lati - linea, 0); nval += ok
        out["terra"][k] = dterra; out["h15"][k] = h15; out["ang"][k] = ang; out["hc"][k] = hc
        out["valle"][k] = val / np.maximum(nval, 1)
    return out


def tabella(nome, x, rapp, err, limiti, unita=""):
    righe = [f"\n{nome}"]
    for lo, hi in zip(limiti[:-1], limiti[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() < 300:
            continue
        righe.append(f"  {lo:>6}-{hi:<6}{unita:3} celle {m.sum():7d}   radar/modello {np.median(rapp[m]):.2f}   errore {np.mean(err[m]):+5.2f} kn")
    return righe


def main():
    cart, ril_p, usc = sys.argv[1:4]
    os.makedirs(usc, exist_ok=True)
    ril = Rilievo(ril_p)
    tutti = {k: [] for k in ("area", "dist", "terra", "h15", "ang", "hc", "valle", "rapp", "err", "vel", "ora", "mese")}
    for area in AREE:
        z, mu, mv = campi(cart, area)
        lat, lon = z["lat"], z["lon"]
        f = caratteristiche(lat, lon, ril)
        R = z["radar"].astype("f8"); M = np.hypot(mu, mv)
        da = (np.degrees(np.arctan2(-mu, -mv)) + 360) % 360                  # da dove viene il vento
        k = (np.round(da / (360 / NSET)).astype(int)) % NSET
        D = np.broadcast_to(z["dist"], R.shape)
        ok = np.isfinite(R) & np.isfinite(M) & (M > 3.0) & (D <= 10)
        si, sj, sk = np.nonzero(ok)
        kk = k[si, sj, sk]
        ora = ((z["t"] // 3600) % 24)[si]
        mese = np.array([int(np.datetime64(int(t), "s").astype("datetime64[M]").astype(int) % 12 + 1) for t in z["t"]])[si]
        tutti["area"].append(np.full(len(si), AREE.index(area)))
        tutti["dist"].append(D[si, sj, sk])
        for nome in ("terra", "h15", "ang", "hc", "valle"):
            tutti[nome].append(f[nome][kk, sj, sk])
        tutti["rapp"].append(R[si, sj, sk] / M[si, sj, sk])
        tutti["err"].append((M[si, sj, sk] - R[si, sj, sk]) * KN)
        tutti["vel"].append(M[si, sj, sk] * KN)
        tutti["ora"].append(ora); tutti["mese"].append(mese)
        print(f"{area}: {len(si)} celle entro 10 km con il modello oltre 6 nodi")
    T = {k: np.concatenate(v) for k, v in tutti.items()}
    np.savez_compressed(os.path.join(usc, "diagnosi_celle.npz"), **T)
    righe = [f"Celle entro 10 km dalla terra, modello letto dal mare oltre 6 nodi: {len(T['rapp'])}",
             f"Tutte: radar/modello {np.median(T['rapp']):.2f}, errore {np.mean(T['err']):+.2f} kn"]
    da_terra = T["terra"] < RMAX
    righe += [f"Vento che arriva da terra (entro 30 km sopravvento): {da_terra.mean() * 100:.0f}% delle celle",
              f"  da terra: radar/modello {np.median(T['rapp'][da_terra]):.2f}, errore {np.mean(T['err'][da_terra]):+.2f} kn",
              f"  dal mare: radar/modello {np.median(T['rapp'][~da_terra]):.2f}, errore {np.mean(T['err'][~da_terra]):+.2f} kn"]
    m = da_terra
    righe += tabella("Vento da terra, per quota massima sopravvento entro 15 km dalla costa (m)", T["h15"][m], T["rapp"][m], T["err"][m], [0, 100, 300, 600, 900, 1300, 3000])
    righe += tabella("Vento da terra, per angolo sotto cui si vede il rilievo (gradi)", T["ang"][m], T["rapp"][m], T["err"][m], [0, 1, 2, 4, 7, 12, 90])
    righe += tabella("Vento da terra, per distanza dalla prima terra sopravvento (km)", T["terra"][m], T["rapp"][m], T["err"][m], [0, 1, 2, 4, 8, 15, 30])
    righe += tabella("Vento da terra, valle o passo sopravvento (m: rilievo ai lati sopra la linea del vento)", T["valle"][m], T["rapp"][m], T["err"][m], [-2000, -100, 0, 100, 250, 500, 3000])
    righe += tabella("Vento da terra, per forza del modello (nodi)", T["vel"][m], T["rapp"][m], T["err"][m], [6, 10, 14, 18, 25, 60])
    righe += tabella("Vento dal mare, per forza del modello (nodi)", T["vel"][~m], T["rapp"][~m], T["err"][~m], [6, 10, 14, 18, 25, 60])
    righe.append("\nPer area (vento da terra | dal mare): radar/modello, errore in nodi")
    for a, nome in enumerate(AREE):
        x = T["area"] == a
        righe.append(f"  {nome:<20} {np.median(T['rapp'][x & m]):.2f} {np.mean(T['err'][x & m]):+5.2f} | {np.median(T['rapp'][x & ~m]):.2f} {np.mean(T['err'][x & ~m]):+5.2f}   ({(x & m).sum()} | {(x & ~m).sum()} celle)")
    righe.append("\nPer ora del passaggio (UTC) e stagione, vento da terra: radar/modello")
    for o in (5, 16, 17):
        x = m & (T["ora"] == o)
        if x.sum() > 300:
            righe.append(f"  ore {o:02d}: {np.median(T['rapp'][x]):.2f} ({x.sum()} celle)")
    est = np.isin(T["mese"], [6, 7, 8, 9]); inv = np.isin(T["mese"], [11, 12, 1, 2, 3])
    righe.append(f"  giu-set: {np.median(T['rapp'][m & est]):.2f}   nov-mar: {np.median(T['rapp'][m & inv]):.2f}")
    testo = "\n".join(righe)
    print(testo)
    open(os.path.join(usc, "diagnosi_rilievo.txt"), "w").write(testo + "\n")


if __name__ == "__main__":
    main()
