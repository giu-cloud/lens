"""Lens 2 sopra il primo Lens: quanto migliora contro il radar e quanto disegno di Lens lascia (26/11).

La domanda di Giano: un solo passaggio che si veda come la terza colonna (primo Lens a 100 m), e altre aggiunte solo se
migliorano senza peggiorare. Qui, sulle settimane della cassaforte (libreria di Lens 2 costruita senza di esse):

  modello      il modello letto dal mare (colonna "mod" della tabella)
  0.54         modello x Lens 0.54
  primo        modello x primo Lens (la terza colonna): Fr 0,7 con aria stabile, Fr 2 con aria mescolata, a pesi
  lens2        modello x fattore della libreria (= Lens 2 relativo cella per cella sopra il primo Lens, che sulla cella
               del radar dà proprio il fattore: di Lens resta solo il disegno dentro la cella)
  scala s      il radar fissa il livello oltre s km, Lens il disegno sotto: modello x primo x G_s / R_s, con G_s e R_s
               il fattore e il primo Lens lisciati a s km (gaussiana con sigma s/2, solo mare)
  peso w       modello x primo x (fattore / primo)^w

Misure (nodi), solo dove il radar vede: errore quadratico, errore sistematico per cella (tutto ed entro 3 km),
errore medio. "Disegno" = quanta parte del disegno di Lens sopra il km si vede ancora: pendenza di log(previsione /
modello) su log(primo Lens), per passaggio, centrate, mediana sui passaggi (1 = tutto, come la terza colonna).

  python3 lens2_su_primo.py <area> <atlante .lna> <campi_area.npz> <libreria senza cassaforte .lra>
                            <primo Fr 0,7 .bin> <primo Fr 2 .bin> <uscita .json>
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

qui = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [qui, os.path.join(qui, "..", "cfd")]
from controlla_atlante import leggi as leggi_lna      # noqa: E402
from leggi_campi import leggi                          # noqa: E402
from libreria_radar import Libreria                    # noqa: E402
from stazioni_cfd import vicini                        # noqa: E402
from taratura_radar import su_celle_radar              # noqa: E402

AREE14 = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova",
          "pontine", "corsica-nord", "monaco", "cicladi"]
SETTIMANA = 7 * 86400
SCALE = tuple(float(x) for x in os.environ.get("SCALE", "3,5,10,20").split(","))
PESI = tuple(float(x) for x in os.environ.get("PESI", "0.25,0.5,0.75").split(",") if x)
DA = np.arange(16) * 22.5


def liscia(Z, mask, sig):
    a0 = np.where(mask, Z, 0.0)
    out = gaussian_filter(a0, sig) / np.maximum(gaussian_filter(mask.astype(float), sig), 1e-6)
    return np.where(mask, out, np.nan)


def risposta(p_bin, mare, rlat, rlon, dist_c):
    la, lo, F = leggi(p_bin, 16)
    S = np.stack([np.hypot(*F[k][1]) for k in range(16)])
    Sr, _ = su_celle_radar(S, mare, la, lo, rlat, rlon)
    Soff = []
    for k in range(16):
        v = np.nan
        for soglia in (8, 6, 4):
            sel = (dist_c > soglia) & np.isfinite(Sr[k])
            if sel.sum() >= 5:
                v = np.nanmedian(Sr[k][sel]); break
        Soff.append(v)
    return Sr / np.array(Soff)[:, None, None]


def main():
    area, p_lna, p_campi, p_lib, p07, p20, uscita = sys.argv[1:8]
    A = leggi_lna(p_lna); m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
    C = np.load(p_campi); rlat, rlon, dist_c = C["lat"].astype(float), C["lon"].astype(float), C["dist"]
    dkm = abs(rlat[1] - rlat[0]) * 111.19
    R07 = risposta(p07, mare, rlat, rlon, dist_c); R20 = risposta(p20, mare, rlat, rlon, dist_c)
    valido = np.isfinite(R07).all(0) & np.isfinite(R20).all(0)
    Rs = {s: (np.stack([liscia(R07[k], valido, s / 2 / dkm) for k in range(16)]),
              np.stack([liscia(R20[k], valido, s / 2 / dkm) for k in range(16)])) for s in SCALE}

    T = np.load(os.environ.get("TABELLA", "/home/claude/lavoro/tabella6.npz"), mmap_mode="r")
    sel = np.asarray(T["area"]) == AREE14.index(area)
    g = lambda k: np.asarray(T[k][sel])                                   # noqa: E731
    j, i = g("j").astype(int), g("i").astype(int)
    vero = np.load(os.environ.get("VERO", "/home/claude/lavoro/radar_vero_tabella6.npy"), mmap_mode="r")[sel]
    mod, D0, S0, dist, t, sc = g("mod"), g("D0"), g("S0"), g("dist"), g("t"), g("scena")
    qL, dT925, dirc = g("q_lens054"), g("dT925"), g("dir")
    v925, d925, ora = g("v925"), g("d925"), g("ora")
    cass = ((t // SETTIMANA) % 5 == 2)
    neu = np.clip((4 - dT925) / 4, 0, 1)                                  # 1 = aria mescolata, 0 = stabile
    neu = np.where(np.isfinite(neu), neu, 0.5)
    k0, k1, w = vicini(DA, D0)

    def su_righe(R):
        return (1 - w) * R[k0, j, i] + w * R[k1, j, i]
    rc = (1 - neu) * su_righe(R07) + neu * su_righe(R20)
    rcs = {s: (1 - neu) * su_righe(Rs[s][0]) + neu * su_righe(Rs[s][1]) for s in SCALE}

    L = Libreria(p_lib)
    assert (L.ny, L.nx) == (len(rlat), len(rlon)), "la griglia della libreria non è quella del radar"
    cella = j * L.nx + i; pos = {c: k for k, c in enumerate(L.celle)}
    jl, il = np.array(L.celle) // L.nx, np.array(L.celle) % L.nx
    mask_lib = np.zeros((L.ny, L.nx), bool); mask_lib[jl, il] = True
    G = np.full(len(j), np.nan); Gs = {s: np.full(len(j), np.nan) for s in SCALE}
    for s_id in np.unique(sc[cass]):
        r = np.where(sc == s_id)[0]
        a = np.radians(D0[r[0]]); b = np.radians(d925[r[0]])
        if not (np.isfinite(a) and np.isfinite(b) and np.isfinite(S0[r[0]]) and np.isfinite(v925[r[0]])):
            continue
        dc = np.full(L.nc, D0[r[0]]); idx = np.array([pos.get(c, -1) for c in cella[r]])
        dc[idx[idx >= 0]] = dirc[r[idx >= 0]]
        gg = L.fattore(S0[r[0]] * np.sin(a), S0[r[0]] * np.cos(a), v925[r[0]] * np.sin(b), v925[r[0]] * np.cos(b),
                       1.0 if ora[r[0]] >= 12 else 0.0, dc)
        ok_i = idx >= 0
        G[r[ok_i]] = gg[idx[ok_i]]
        griglia = np.full((L.ny, L.nx), np.nan); griglia[jl, il] = np.log(gg)
        for s in SCALE:
            ls = liscia(griglia, mask_lib & np.isfinite(griglia), s / 2 / dkm)
            Gs[s][r] = np.exp(ls[j[r], i[r]])

    use = (cass & np.isfinite(vero) & (vero > 0.5) & (mod > 2.0) & np.isfinite(rc) & (rc > 0) & np.isfinite(G)
           & np.isfinite(qL) & (qL > 0) & (S0 > 2.0) & np.all([np.isfinite(rcs[s]) & np.isfinite(Gs[s]) for s in SCALE], 0))
    cel = j * 1000 + i

    def sist(e, mk):
        ids, inv = np.unique(cel[mk], return_inverse=True); n = np.bincount(inv); me = np.bincount(inv, e[mk]) / n
        return float(np.sqrt(np.sum(n * me ** 2) / n.sum()))

    def disegno(pred):
        """Pendenza di log(pred/mod) su log(primo), centrate per passaggio; mediana sui passaggi con almeno 30 celle."""
        x = np.log(rc); y = np.log(pred / mod); out = []
        for s_id in np.unique(sc[use]):
            k = use & (sc == s_id)
            if k.sum() < 30:
                continue
            xc = x[k] - x[k].mean(); yc = y[k] - y[k].mean()
            if np.sum(xc ** 2) > 1e-6:
                out.append(np.sum(xc * yc) / np.sum(xc ** 2))
        return float(np.median(out)) if out else None

    varianti = {"modello": mod, "0.54": mod * qL, "primo (terza colonna)": mod * rc, "lens2 (relativo per cella)": mod * G}
    for s in SCALE:
        varianti[f"scala {s:g} km"] = mod * rc * Gs[s] / rcs[s]
    for p in PESI:
        varianti[f"peso {p:g}"] = mod * rc * (G / rc) ** p
    vic = dist < 3
    ris = {"area": area, "righe": int(use.sum()), "passaggi": int(len(np.unique(sc[use]))), "varianti": {}}
    print(f"{area}: {use.sum()} righe della cassaforte, {len(np.unique(sc[use]))} passaggi")
    for nome, pred in varianti.items():
        e = (pred - vero) / 0.5144
        d = {"rmse": float(np.sqrt(np.mean(e[use] ** 2))), "sist": sist(e, use), "media": float(np.mean(e[use])),
             "sist_3km": sist(e, use & vic), "media_3km": float(np.mean(e[use & vic])) if (use & vic).any() else None,
             "disegno": disegno(pred)}
        ris["varianti"][nome] = d
        print(f"  {nome:28s} rmse {d['rmse']:.2f}  sist {d['sist']:.2f}  media {d['media']:+.2f}  | <3 km sist "
              f"{d['sist_3km']:.2f} media {d['media_3km']:+.2f}  | disegno {d['disegno']:.2f}", flush=True)
    json.dump(ris, open(uscita, "w"), indent=1)


if __name__ == "__main__":
    main()
