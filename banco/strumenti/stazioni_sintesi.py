"""Riassunto della prova con le stazioni: errore di ECMWF e ICON-2I, grezzi e corretti, e del loro insieme, contro il vento
sul mare davanti alle stazioni buone (stazione divisa per la sua esposizione nel settore), nelle settimane della cassaforte.

  python3 stazioni_sintesi.py <cartella con i json di stazioni_prova> <uscita.json>

Due modi di giudicare:
- assoluto: errore medio e quadratico contro la stazione riportata al mare. Prova le correzioni della forza;
- libero: per ogni stazione e previsione la scala migliore (un solo numero per tutta la serie) prima dell'errore
  quadratico. Toglie qualunque errore costante e misura solo se la correzione segue meglio i cambiamenti.
Le ore: 00-03 notte, 06-09 mattina, 12-15 pomeriggio (brezze, lontano dai passaggi del radar), 18-21 sera.
"""
import glob
import json
import os
import sys

import numpy as np

KN = 1.943844
FASCE = {"notte (00-03)": (0, 3), "mattina (06-09)": (6, 9), "pomeriggio (12-15)": (12, 15), "sera (18-21)": (18, 21)}
P_E, P_I = 0.6, 0.4


def metriche(oss, prev, libero=False):
    ok = np.isfinite(oss) & np.isfinite(prev)
    o, p = oss[ok], prev[ok]
    if len(o) < 20:
        return None
    if libero:
        p = p * np.sum(o * p) / np.sum(p * p)
    return dict(n=int(len(o)), medio=round(float(np.mean(p - o)) * KN, 2), quadratico=round(float(np.sqrt(np.mean((p - o) ** 2))) * KN, 2))


def tabella(R, libero=False):
    o = R["staz_mare"]
    var = {"ECMWF": R["ecmwf"], "ECMWF corretto": R["ecmwf_c"]}
    if np.isfinite(R["icon"]).any():
        var.update({"ICON-2I": R["icon"], "ICON-2I corretto": R["icon_c"],
                    "insieme grezzo": P_E * R["ecmwf"] + P_I * R["icon"], "insieme corretto": P_E * R["ecmwf_c"] + P_I * R["icon_c"]})
    return {k: metriche(o, v, libero) for k, v in var.items()}


def libero_per_stazione(R, maschera):
    """Errore quadratico dopo la scala migliore per ogni stazione: media pesata sulle stazioni."""
    out = {}
    for nome in ("ecmwf", "ecmwf_c", "icon", "icon_c", "ens", "ens_c"):
        num = den = 0.0
        for st in np.unique(R["st"][maschera]):
            m = maschera & (R["st"] == st)
            o = R["staz_mare"][m]; p = R[nome][m]
            ok = np.isfinite(o) & np.isfinite(p)
            if ok.sum() < 20:
                continue
            o, p = o[ok], p[ok]
            p = p * np.sum(o * p) / np.sum(p * p)
            num += np.sum((p - o) ** 2); den += len(o)
        out[nome] = round(float(np.sqrt(num / den)) * KN, 3) if den else None
    return out


def main():
    cart, usc = sys.argv[1:3]
    righe, stazioni = [], []
    for f in sorted(glob.glob(os.path.join(cart, "*.json"))):
        d = json.load(open(f))
        if "righe" not in d:
            continue
        for r in d["righe"]:
            r["area"] = d["area"]
        righe += d["righe"]
        stazioni += [dict(s, area=d["area"]) for s in d["stazioni"]]
    chiavi = ["staz", "staz_mare", "ecmwf", "ecmwf_c", "icon", "icon_c", "disp", "S0"]
    R = {k: np.array([r.get(k, np.nan) for r in righe], float) for k in chiavi}
    R["st"] = np.array([r["st"] for r in righe]); R["area"] = np.array([r["area"] for r in righe])
    R["ora"] = np.array([r["ora"] for r in righe])
    R["ens"] = P_E * R["ecmwf"] + P_I * R["icon"]; R["ens_c"] = P_E * R["ecmwf_c"] + P_I * R["icon_c"]
    sopra = R["S0"] >= 2.0                       # con il vento al largo sotto 4 nodi i rapporti sono solo rumore
    tutte = sopra & np.isfinite(R["icon"])        # stesse ore per tutti i modelli
    out = {"stazioni_usate": sorted(set(R["st"].tolist())), "confronti": int(sopra.sum()),
           "con_icon": int(tutte.sum()), "assoluto": {}, "libero": {}, "per_stazione": {}, "per_ora": {}}
    sub = {k: (v[tutte] if isinstance(v, np.ndarray) else v) for k, v in R.items()}
    out["assoluto"]["tutte le ore"] = tabella(sub)
    out["libero"]["tutte le ore"] = libero_per_stazione(R, tutte)
    for nome, (h0, h1) in FASCE.items():
        m = tutte & (R["ora"] >= h0) & (R["ora"] <= h1)
        out["assoluto"][nome] = tabella({k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in R.items()})
        out["libero"][nome] = libero_per_stazione(R, m)
    for st in sorted(set(R["st"].tolist())):
        m = sopra & (R["st"] == st)
        s = next((x for x in stazioni if f"{x['fonte']}_{x['codice']}" == st), {})
        out["per_stazione"][st] = dict(nome=s.get("nome"), area=s.get("area"), ottanti=s.get("ottanti_buoni"),
                                       fattore_medio=round(float(np.nanmean(R["ecmwf_c"][m] / R["ecmwf"][m])), 3),
                                       assoluto=tabella({k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in R.items()}),
                                       libero=libero_per_stazione(R, m))
    # fiducia: dispersione degli analoghi contro errore relativo dell'ECMWF corretto
    m = sopra & np.isfinite(R["disp"]) & (R["staz_mare"] > 1)
    if m.sum() > 100:
        q = np.quantile(R["disp"][m], [0.2, 0.4, 0.6, 0.8])
        cl = np.digitize(R["disp"][m], q)
        rel = np.abs(R["ecmwf_c"][m] - R["staz_mare"][m]) / np.maximum(R["staz_mare"][m], 2.0)
        out["fiducia"] = [dict(fino_a=round(float(q[i]), 3) if i < 4 else None, errore_relativo=round(float(np.median(rel[cl == i])), 3),
                               n=int((cl == i).sum())) for i in range(5)]
    out["scartate"] = {f"{s['fonte']}_{s['codice']}": dict(nome=s["nome"], area=s["area"], motivo=s["scartata"])
                       for s in stazioni if s.get("scartata")}
    out["tenute"] = {f"{s['fonte']}_{s['codice']}": dict(nome=s["nome"], area=s["area"], celle=s["celle"],
                                                          ottanti=s["ottanti_buoni"], taratura=s["taratura"])
                     for s in stazioni if not s.get("scartata")}
    json.dump(out, open(usc, "w"), indent=1, ensure_ascii=False)
    print(f"{len(out['stazioni_usate'])} stazioni, {out['confronti']} confronti ({out['con_icon']} con ICON-2I)")
    for chiave in ("assoluto",):
        print(f"\n== {chiave} (nodi: errore medio / quadratico)")
        for fascia, t in out[chiave].items():
            print(f"  {fascia:20s} " + "  ".join(f"{k}: {v['medio']:+.2f}/{v['quadratico']:.2f}" for k, v in t.items() if v))
    print("\n== libero (errore quadratico dopo la scala migliore per stazione, nodi)")
    for fascia, t in out["libero"].items():
        print(f"  {fascia:20s} " + "  ".join(f"{k}: {v}" for k, v in t.items()))
    print("\n== per stazione (assoluto: ECMWF -> corretto; insieme grezzo -> corretto)")
    for st, v in out["per_stazione"].items():
        a = v["assoluto"]
        e, ec = a.get("ECMWF"), a.get("ECMWF corretto")
        eg, ecc = a.get("insieme grezzo"), a.get("insieme corretto")
        print(f"  {st:22s} {str(v['nome'])[:16]:16s} {v['area']:20s} fattore {v['fattore_medio']:.2f}  "
              f"ECMWF {e['medio']:+.2f}/{e['quadratico']:.2f} -> {ec['medio']:+.2f}/{ec['quadratico']:.2f}"
              + (f"   insieme {eg['medio']:+.2f}/{eg['quadratico']:.2f} -> {ecc['medio']:+.2f}/{ecc['quadratico']:.2f}" if eg and ecc else "")
              + f"   n={e['n']}")
    if "fiducia" in out:
        print("\n== fiducia", out["fiducia"])


if __name__ == "__main__":
    main()
