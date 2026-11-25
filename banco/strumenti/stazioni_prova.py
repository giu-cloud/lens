"""Prova di Lens 2 contro le stazioni a terra, a tutte le ore, con le stazioni scelte e tarate dal radar.

  python3 stazioni_prova.py <area> <elenco stazioni .json> <cartella stazioni> <librerie senza cassaforte> \
      <ecmwf passi 3-6> <ecmwf passi 9-12> <stab 3-6> <925 9-12> <icon2i ogni 3 ore .npz> <rilievo> <uscita.json>

Il problema delle stazioni a terra è che leggono il vento del posto (città, alberi, collina), non quello del mare.
Qui il radar fa da giudice, settore per settore:
1. Per ogni stazione si prendono le celle della libreria entro R km (il mare davanti alla stazione).
2. Ai passaggi del radar delle settimane di addestramento (le stesse delle librerie) si confronta la stazione con il
   radar medio su quelle celle, per ottante della direzione del vento al largo. Un ottante è "buono" se la stazione va
   d'accordo con il mare: almeno N passaggi, correlazione alta, rapporto stazione/mare stabile e non lontano da 1.
   Il rapporto medio di ogni ottante buono è l'esposizione della stazione in quel settore.
3. Nelle settimane della cassaforte, a tutte le ore (ogni 3 ore), la stazione divisa per la sua esposizione diventa una
   misura del vento sul mare davanti, e con quella si giudicano ECMWF e ICON-2I, grezzi e corretti, e il loro insieme.
   Gli ottanti non buoni si scartano: una stazione di città vale solo per i venti con cui vede il mare.
"""
import datetime as dt
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from banco_icon2i import terra_nodi
from boe import lettura_ecmwf
from libreria_radar import Libreria
from rilievo import Rilievo
from stabilita import leggi

SETTIMANA = 7 * 86400
R_KM = float(os.environ.get("R_KM", "5"))
MIN_CELLE = 3
MIN_N, MIN_COR, RAPP, MAX_SCARTO = 12, 0.6, (0.5, 1.6), 0.40
VMIN = 2.5          # m/s sul radar per tarare
QUOTA_MAX = 200.0   # m: sopra, la stazione sente un altro vento (capi alti, colli)
FINESTRA = 20 * 60  # s tra misura e ora del modello


def cassa(t):
    return (np.asarray(t) // SETTIMANA) % 5 == 2


def ottante(d):
    return int(np.round(d / 45.0)) % 8


def km(lat1, lon1, lat2, lon2):
    return 111.2 * np.hypot(lat1 - lat2, (lon1 - lon2) * np.cos(np.radians(lat2)))


def celle_vicine(L, lat, lon):
    """Celle della libreria entro R_KM dalla stazione, lasciando fuori la fascia di 1,5 km dalla costa, dove il radar
    è spesso vuoto o sporcato dalla terra."""
    la = L.lat0 + (L.celle // L.nx) * L.dlat; lo = L.lon0 + (L.celle % L.nx) * L.dlon
    d = km(la, lo, lat, lon)
    return np.nonzero((d <= R_KM) & (L.dist >= 1.5))[0], la, lo


def misura_vicina(S, t, finestra=FINESTRA):
    k = np.searchsorted(S["t"], t)
    best = None
    for j in (k - 1, k):
        if 0 <= j < len(S["t"]) and abs(S["t"][j] - t) <= finestra and np.isfinite(S["vento"][j]):
            if best is None or abs(S["t"][j] - t) < abs(S["t"][best] - t):
                best = j
    return best


def taratura(S, LE, cel):
    """Esposizione per ottante dai passaggi del radar nella libreria (settimane di addestramento)."""
    coppie = []
    for k in range(LE.ns):
        R = LE.M[k, cel] * LE.Q[k, cel]
        if np.isfinite(R).sum() < max(1, len(cel) // 3):
            continue
        j = misura_vicina(S, LE.t[k], 30 * 60)
        if j is None:
            continue
        D0 = (np.degrees(np.arctan2(LE.a0[k], LE.b0[k])) + 360) % 360
        coppie.append((ottante(D0), float(S["vento"][j]), float(np.nanmean(R))))
    tab = {}
    for o in range(8):
        x = np.array([(s, r) for oo, s, r in coppie if oo == o and r >= VMIN])
        if len(x) < 3:
            tab[o] = dict(n=int(len(x)), buono=False); continue
        s, r = x[:, 0], x[:, 1]
        a = float(np.sum(s * r) / np.sum(r * r))
        cor = float(np.corrcoef(s, r)[0, 1]) if np.std(s) > 0 and np.std(r) > 0 else 0.0
        scarto = float(np.sqrt(np.mean((s - a * r) ** 2)) / np.mean(r))
        buono = len(x) >= MIN_N and cor >= MIN_COR and RAPP[0] <= a <= RAPP[1] and scarto <= MAX_SCARTO
        tab[o] = dict(n=int(len(x)), rapporto=round(a, 3), cor=round(cor, 3), scarto=round(scarto, 3), buono=bool(buono))
    return tab, len(coppie)


def main():
    area, elenco_p, cart_s, cart_l, e36, e912, s36, p912, icon_p, ril_p, usc = sys.argv[1:12]
    elenco = [s for s in json.load(open(elenco_p)) if s["area"] == area]
    LE = Libreria(os.path.join(cart_l, f"lra1_ecmwf_{area}.lra"))
    LI_p = os.path.join(cart_l, f"lra1_icon-2i_{area}.lra")
    LI = Libreria(LI_p) if os.path.exists(LI_p) else None
    ril = Rilievo(ril_p)
    stazioni = []
    for s in elenco:
        f = os.path.join(cart_s, f"{s['fonte']}_{s['codice']}.npz")
        if not os.path.exists(f):
            continue
        cel, la, lo = celle_vicine(LE, s["lat"], s["lon"])
        info = dict(s, celle=int(len(cel)), dist_mare_km=round(float(km(la, lo, s["lat"], s["lon"]).min()), 2))
        if np.isfinite(s.get("quota", np.nan)) and s["quota"] > QUOTA_MAX:
            info["scartata"] = f"quota {s['quota']:.0f} m"; stazioni.append((info, None, None, None)); continue
        if len(cel) < MIN_CELLE:
            info["scartata"] = f"poche celle di mare entro {R_KM:.0f} km"; stazioni.append((info, None, None, None)); continue
        Z = np.load(f); o = np.argsort(Z["t"])
        S = {k: Z[k][o] for k in ("t", "vento", "direzione")}
        tab, n = taratura(S, LE, cel)
        info["taratura"] = tab; info["passaggi"] = n
        buoni = [o_ for o_, v in tab.items() if v["buono"]]
        info["ottanti_buoni"] = buoni
        if not buoni:
            info["scartata"] = "nessun settore d'accordo con il radar"; stazioni.append((info, None, None, None)); continue
        celI = celle_vicine(LI, s["lat"], s["lon"])[0] if LI is not None else None
        stazioni.append((info, S, cel, celI))
    attive = [x for x in stazioni if x[1] is not None]
    print(f"{area}: {len(elenco)} stazioni, {len(attive)} con almeno un settore buono", flush=True)
    for info, *_ in stazioni:
        print(f"  {info['fonte']} {info['codice']} {info['nome'][:28]:28s} celle {info['celle']:3d} "
              f"{info.get('scartata') or 'ottanti ' + ','.join(str(o) for o in info['ottanti_buoni'])}", flush=True)
    righe = []
    if attive:
        lat_c = LE.lat0 + (LE.celle // LE.nx) * LE.dlat; lon_c = LE.lon0 + (LE.celle % LE.nx) * LE.dlon
        if LI is not None and os.path.exists(icon_p):
            Zi = np.load(icon_p); it = Zi["t"]; oi = np.argsort(it); it = it[oi]
            IU, IV = Zi["u"][oi], Zi["v"][oi]
            pesi_i = lettura_ecmwf(Zi["lat"], Zi["lon"], terra_nodi(ril, Zi["lat"], Zi["lon"]),
                                   LI.lat0 + (LI.celle // LI.nx) * LI.dlat, LI.lon0 + (LI.celle % LI.nx) * LI.dlon)
        else:
            it = None
        pesi = None
        tutte_t = np.unique(np.concatenate([np.round(x[1]["t"] / 3600.0) * 3600 for x in attive]).astype(np.int64))
        tutte_t = tutte_t[(tutte_t % 10800 == 0) & cassa(tutte_t)]
        cache = {}
        for tv in tutte_t:
            d = dt.datetime.fromtimestamp(int(tv), dt.timezone.utc)
            run = d.replace(hour=0 if 3 <= d.hour <= 12 else 12, minute=0, second=0)
            if d.hour == 0:
                run = run - dt.timedelta(days=1)
            passo = int((d - run).total_seconds() // 3600)
            nome = f"ecmwf_{run:%Y%m%d%H}.grb2"
            f_w = os.path.join(e36, nome) if passo <= 6 else os.path.join(e912, nome.replace("ecmwf_", "e912_"))
            f_s = os.path.join(s36, nome.replace("ecmwf_", "stab_")) if passo <= 6 else os.path.join(p912, nome.replace("ecmwf_", "p912_"))
            if not os.path.exists(f_w):
                continue
            if f_w not in cache:
                cache.clear(); cache[f_w] = leggi(f_w)
            lat, lon, C = cache[f_w]
            if ("10u", passo) not in C:
                continue
            if pesi is None:
                LA, LO = np.meshgrid(lat, lon, indexing="ij")
                terra = np.zeros(LA.shape, bool)
                c_la, c_lo = lat_c.mean(), lon_c.mean()
                for j in range(LA.shape[0]):
                    for i in range(LA.shape[1]):
                        if abs(LA[j, i] - c_la) < 1.5 and abs(LO[j, i] - c_lo) < 2:
                            _, _, hh = ril.griglia(0, LA[j, i] - 0.003, LA[j, i] + 0.003, LO[j, i] - 0.004, LO[j, i] + 0.004)
                            terra[j, i] = hh.max() > 0
                pesi = lettura_ecmwf(lat, lon, terra, lat_c, lon_c)
                nodi_area = (LA >= lat_c.min() - 0.25) & (LA <= lat_c.max() + 0.25) & (LO >= lon_c.min() - 0.25) & (LO <= lon_c.max() + 0.25)
            U = C[("10u", passo)]; V = C[("10v", passo)]
            u = np.array([sum(w * U[j, i] for j, i, w in p) for p in pesi], "f4"); v = np.array([sum(w * V[j, i] for j, i, w in p) for p in pesi], "f4")
            a925 = b925 = None
            if os.path.exists(f_s):
                _, _, Sx = leggi(f_s)
                if ("u@925", passo) in Sx:
                    u9 = float(Sx[("u@925", passo)][nodi_area].mean()); v9 = float(Sx[("v@925", passo)][nodi_area].mean())
                    s9 = np.hypot(u9, v9); d9 = np.arctan2(-u9, -v9); a925, b925 = s9 * np.sin(d9), s9 * np.cos(d9)
            largo = LE.dist >= 8
            if largo.sum() < 20:
                largo = LE.dist >= np.quantile(LE.dist[LE.dist > 0], 0.75)
            S0 = float(np.median(np.hypot(u, v)[largo])); D0r = np.arctan2(-u[largo].mean(), -v[largo].mean())
            a0, b0 = S0 * np.sin(D0r), S0 * np.cos(D0r)
            if a925 is None:
                a925, b925 = a0, b0
            sera = 1.0 if 12 <= d.hour <= 23 else 0.0
            dirc = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
            G, SP = LE.fattore(a0, b0, a925, b925, sera, dirc, dispersione=True)
            oc = ottante((np.degrees(D0r) + 360) % 360)
            vel = np.hypot(u, v)
            gi = None
            if it is not None:
                kk = np.searchsorted(it, tv)
                if kk < len(it) and it[kk] == tv:
                    Uk = IU[kk].astype("f4"); Vk = IV[kk].astype("f4")
                    ui = np.array([sum(w * Uk[j, i] for j, i, w in pp) for pp in pesi_i], "f4")
                    vi = np.array([sum(w * Vk[j, i] for j, i, w in pp) for pp in pesi_i], "f4")
                    li = LI.dist >= 8
                    if li.sum() < 20:
                        li = LI.dist >= np.quantile(LI.dist[LI.dist > 0], 0.75)
                    S0i = float(np.median(np.hypot(ui, vi)[li])); D0i = np.arctan2(-ui[li].mean(), -vi[li].mean())
                    dici = (np.degrees(np.arctan2(-ui, -vi)) + 360) % 360
                    gi = (np.hypot(ui, vi), LI.fattore(S0i * np.sin(D0i), S0i * np.cos(D0i), a925, b925, sera, dici))
            for info, S, cel, celI in attive:
                tb = info["taratura"][oc]
                if not tb["buono"]:
                    continue
                j = misura_vicina(S, tv)
                if j is None:
                    continue
                r = {"st": info["fonte"] + "_" + info["codice"], "t": int(tv), "ora": d.hour, "ott": oc, "S0": round(S0, 2),
                     "staz": float(S["vento"][j]), "staz_mare": float(S["vento"][j] / tb["rapporto"]),
                     "ecmwf": float(np.mean(vel[cel])), "ecmwf_c": float(np.mean(vel[cel] * G[cel])),
                     "disp": float(np.nanmean(SP[cel])), "passo": passo}
                if gi is not None and celI is not None and len(celI):
                    r["icon"] = float(np.mean(gi[0][celI])); r["icon_c"] = float(np.mean(gi[0][celI] * gi[1][celI]))
                righe.append(r)
    out = {"area": area, "stazioni": [x[0] for x in stazioni], "righe": righe,
           "regole": dict(R_KM=R_KM, MIN_CELLE=MIN_CELLE, MIN_N=MIN_N, MIN_COR=MIN_COR, RAPP=RAPP, MAX_SCARTO=MAX_SCARTO,
                          VMIN=VMIN, QUOTA_MAX=QUOTA_MAX)}
    json.dump(out, open(usc, "w"))
    print(f"{area}: {len(righe)} confronti -> {usc}", flush=True)


if __name__ == "__main__":
    main()
