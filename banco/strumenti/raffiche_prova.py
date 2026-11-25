"""La raffica sotto Lens 2: va moltiplicata per il fattore della forza o no? Prova con le raffiche delle stazioni.

  python3 raffiche_prova.py <area> <elenco stazioni .json> <cartella stazioni> <librerie senza cassaforte> \
      <ecmwf passi 3-6> <stab 3-6> <rilievo> <uscita.json>

Il radar vede il vento medio, non la raffica: il fattore di Lens 2 è imparato sulla forza media. Qui si guarda se la
raffica misurata segue il fattore come la forza media.
- Raffica del modello: 10fg di ECMWF (massima delle ultime 3 ore), passi 3 e 6 (valide alle 03, 06, 15 e 18 UTC),
  letta dai soli punti di mare sulle celle della stazione, come la forza.
- Raffica della stazione: la massima delle misure nelle 3 ore prima (Météo-France: raffica istantanea massima di ogni
  ora a 10 m). Forza della stazione: quella all'ora del modello.
- Solo le stazioni e gli ottanti buoni della prova con il radar (stazioni_prova.py), solo le settimane della cassaforte.
"""
import datetime as dt
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boe import lettura_ecmwf
from libreria_radar import Libreria
from rilievo import Rilievo
from stabilita import leggi
import stazioni_prova as P


def main():
    area, elenco_p, cart_s, cart_l, e36, s36, ril_p, usc = sys.argv[1:9]
    elenco = [s for s in json.load(open(elenco_p)) if s["area"] == area and s["fonte"] == "mf"]
    LE = Libreria(os.path.join(cart_l, f"lra1_ecmwf_{area}.lra"))
    ril = Rilievo(ril_p)
    attive = []
    for s in elenco:
        f = os.path.join(cart_s, f"{s['fonte']}_{s['codice']}.npz")
        if not os.path.exists(f) or (np.isfinite(s.get("quota", np.nan)) and s["quota"] > P.QUOTA_MAX):
            continue
        cel, _, _ = P.celle_vicine(LE, s["lat"], s["lon"])
        if len(cel) < P.MIN_CELLE:
            continue
        Z = np.load(f); o = np.argsort(Z["t"])
        S = {k: Z[k][o] for k in ("t", "vento", "direzione", "raffica")}
        tab, _ = P.taratura(S, LE, cel)
        if any(v["buono"] for v in tab.values()):
            attive.append((s, S, cel, tab))
    print(f"{area}: {len(attive)} stazioni con raffica e settori buoni: {[s['nome'] for s, *_ in attive]}", flush=True)
    lat_c = LE.lat0 + (LE.celle // LE.nx) * LE.dlat; lon_c = LE.lon0 + (LE.celle % LE.nx) * LE.dlon
    pesi = None; righe = []
    tempi = np.unique(np.concatenate([np.round(x[1]["t"] / 3600.0) * 3600 for x in attive]).astype(np.int64)) if attive else []
    for tv in tempi:
        if not P.cassa(tv):
            continue
        d = dt.datetime.fromtimestamp(int(tv), dt.timezone.utc)
        if d.hour not in (3, 6, 15, 18):
            continue
        run = d.replace(hour=0 if d.hour <= 6 else 12, minute=0, second=0)
        passo = int((d - run).total_seconds() // 3600)
        nome = f"ecmwf_{run:%Y%m%d%H}.grb2"
        f_w = os.path.join(e36, nome); f_s = os.path.join(s36, nome.replace("ecmwf_", "stab_"))
        if not (os.path.exists(f_w) and os.path.exists(f_s)):
            continue
        lat, lon, C = leggi(f_w); _, _, Sx = leggi(f_s)
        if ("10u", passo) not in C or ("10fg", passo) not in Sx:
            continue
        if pesi is None:
            LA, LO = np.meshgrid(lat, lon, indexing="ij")
            terra = np.zeros(LA.shape, bool)
            for j in range(LA.shape[0]):
                for i in range(LA.shape[1]):
                    if abs(LA[j, i] - lat_c.mean()) < 1.5 and abs(LO[j, i] - lon_c.mean()) < 2:
                        _, _, hh = ril.griglia(0, LA[j, i] - 0.003, LA[j, i] + 0.003, LO[j, i] - 0.004, LO[j, i] + 0.004)
                        terra[j, i] = hh.max() > 0
            pesi = lettura_ecmwf(lat, lon, terra, lat_c, lon_c)
            nodi_area = (LA >= lat_c.min() - 0.25) & (LA <= lat_c.max() + 0.25) & (LO >= lon_c.min() - 0.25) & (LO <= lon_c.max() + 0.25)
        U, V, FG = C[("10u", passo)], C[("10v", passo)], Sx[("10fg", passo)]
        u = np.array([sum(w * U[j, i] for j, i, w in p) for p in pesi], "f4"); v = np.array([sum(w * V[j, i] for j, i, w in p) for p in pesi], "f4")
        fg = np.array([sum(w * FG[j, i] for j, i, w in p) for p in pesi], "f4")
        u9 = float(Sx[("u@925", passo)][nodi_area].mean()); v9 = float(Sx[("v@925", passo)][nodi_area].mean())
        s9 = np.hypot(u9, v9); d9 = np.arctan2(-u9, -v9)
        largo = LE.dist >= 8
        if largo.sum() < 20:
            largo = LE.dist >= np.quantile(LE.dist[LE.dist > 0], 0.75)
        S0 = float(np.median(np.hypot(u, v)[largo])); D0r = np.arctan2(-u[largo].mean(), -v[largo].mean())
        G = LE.fattore(S0 * np.sin(D0r), S0 * np.cos(D0r), s9 * np.sin(d9), s9 * np.cos(d9), 1.0 if d.hour >= 12 else 0.0,
                       (np.degrees(np.arctan2(-u, -v)) + 360) % 360)
        oc = P.ottante((np.degrees(D0r) + 360) % 360)
        vel = np.hypot(u, v)
        for s, S, cel, tab in attive:
            if not tab[oc]["buono"]:
                continue
            j = P.misura_vicina(S, tv)
            fin = (S["t"] > tv - 3 * 3600 + 60) & (S["t"] <= tv + 60)
            g3 = S["raffica"][fin]
            if j is None or fin.sum() < 2 or not np.isfinite(g3).any():
                continue
            righe.append({"st": s["codice"], "nome": s["nome"], "ott": oc, "ora": d.hour, "S0": S0,
                          "staz": float(S["vento"][j]), "staz_raff": float(np.nanmax(g3)),
                          "mod": float(np.mean(vel[cel])), "mod_raff": float(np.mean(fg[cel])),
                          "G": float(np.mean(vel[cel] * G[cel]) / max(np.mean(vel[cel]), 0.1)),
                          "G_raff": float(np.mean(fg[cel] * G[cel]) / max(np.mean(fg[cel]), 0.1))})
    json.dump(righe, open(usc, "w"))
    print(f"{area}: {len(righe)} confronti -> {usc}", flush=True)


if __name__ == "__main__":
    main()
