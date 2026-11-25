"""Stazioni a terra con il vento vicino alle aree del banco: aeroporti e stazioni sinottiche (NOAA ISD),
Météo-France (Corsica e Alpi Marittime), più una sonda per le altre fonti.

  python stazioni.py --isd --out stazioni          elenco ISD (i dati ISD sono stati ritirati nel 2025: serve solo l'elenco)
  python stazioni.py --ghcnh --out stazioni        aeroporti e stazioni sinottiche da GHCNh (anche fari e capi dell'Aeronautica)
  python stazioni.py --rmn --out stazioni          anemometri dei porti, Rete Mareografica Nazionale di ISPRA (dati aperti, CC BY 4.0)
  python stazioni.py --francia --out stazioni      stazioni orarie di Météo-France, dipartimenti 20 e 06
  python stazioni.py --sonda --out stazioni        salva le pagine di documentazione delle altre fonti

Uscita: <fonte>_<codice>.npz con t (s UTC), direzione (gradi, NaN se manca o variabile), vento (m/s), raffica (m/s, se c'è),
tipo (codice del messaggio); e stazioni_<fonte>.json con nome, lat, lon, quota, area.
Si tengono le stazioni entro 0,2° dai riquadri delle aree e i dati dal primo aprile 2024.
"""
import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from costa import AREE

MARGINE = 0.2
DA = "2024-04-01"
LOG = []


def log(x):
    LOG.append(str(x)); print(x, flush=True)


def scarica(url, tentativi=4, attesa=5):
    for k in range(tentativi):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "lens-banco/1.0 (ricerca, giu-cloud/lens)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            log(f"  {url}: HTTP {e.code}, riprovo")
        except Exception as e:
            log(f"  {url}: {type(e).__name__} {e}, riprovo")
        time.sleep(attesa * (k + 1))
    return None


def area_di(lat, lon):
    for nome, (s, n, w, e) in AREE.items():
        if s - MARGINE <= lat <= n + MARGINE and w - MARGINE <= lon <= e + MARGINE:
            return nome
    return None


def salva(out, fonte, codice, t, d, v, g, tipo):
    o = np.argsort(t)
    t, d, v, g, tipo = t[o], d[o], v[o], g[o], tipo[o]
    tieni = np.r_[True, np.diff(t) > 0]          # una sola misura per istante
    np.savez_compressed(os.path.join(out, f"{fonte}_{codice}.npz"), t=t[tieni], direzione=d[tieni], vento=v[tieni],
                        raffica=g[tieni], tipo=tipo[tieni])
    return int(tieni.sum())


# ---------------------------------------------------------------- NOAA ISD (aeroporti e sinottiche)

def isd(out):
    testo = scarica("https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv")
    if not testo:
        log("isd-history non raggiungibile"); return
    righe = list(csv.DictReader(io.StringIO(testo.decode("latin-1"))))
    scelte = []
    for r in righe:
        try:
            lat, lon = float(r["LAT"]), float(r["LON"])
        except ValueError:
            continue
        if r["END"] < "20250101":
            continue
        a = area_di(lat, lon)
        if a:
            scelte.append(dict(codice=r["USAF"] + r["WBAN"], nome=r["STATION NAME"].strip(), icao=r["ICAO"].strip(),
                               lat=lat, lon=lon, quota=float(r["ELEV(M)"] or "nan"), area=a,
                               inizio=r["BEGIN"], fine=r["END"]))
    log(f"ISD: {len(scelte)} stazioni nelle aree")
    anni = range(int(DA[:4]), time.gmtime().tm_year + 1)
    t0 = time.mktime(time.strptime(DA, "%Y-%m-%d")) - time.timezone
    for s in scelte:
        T, D, V, G, K = [], [], [], [], []
        for anno in anni:
            b = scarica(f"https://www.ncei.noaa.gov/data/global-hourly/access/{anno}/{s['codice']}.csv")
            if not b:
                continue
            for r in csv.DictReader(io.StringIO(b.decode("latin-1"))):
                w = r.get("WND", "")
                p = w.split(",")
                if len(p) < 5:
                    continue
                dd, qd, ty, ff, qf = p[:5]
                if ff == "9999" or qf in ("2", "3", "6", "7"):
                    continue
                try:
                    t = time.mktime(time.strptime(r["DATE"][:16], "%Y-%m-%dT%H:%M")) - time.timezone
                except ValueError:
                    continue
                if t < t0:
                    continue
                v = int(ff) / 10.0
                d = float(dd) if dd != "999" and qd not in ("2", "3", "6", "7") and ty not in ("V", "9") else np.nan
                if ty == "C":
                    v, d = 0.0, np.nan
                g = np.nan
                gu = r.get("OC1", "")
                if gu:
                    q = gu.split(",")
                    if len(q) >= 2 and q[0] != "9999" and q[1] not in ("2", "3", "6", "7"):
                        g = int(q[0]) / 10.0
                T.append(t); D.append(d); V.append(v); G.append(g); K.append(r.get("REPORT_TYPE", "").strip())
        if not T:
            log(f"  {s['codice']} {s['nome']}: nessun dato"); s["misure"] = 0; continue
        s["misure"] = salva(out, "isd", s["codice"], np.array(T, float), np.array(D, np.float32),
                            np.array(V, np.float32), np.array(G, np.float32), np.array(K))
        log(f"  {s['codice']} {s['nome']} ({s['area']}): {s['misure']} misure")
    json.dump(scelte, open(os.path.join(out, "stazioni_isd.json"), "w"), indent=1, ensure_ascii=False)


# ---------------------------------------------------------------- NOAA GHCNh (ha preso il posto di ISD nel 2025)

GH = "https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly"


def ghcnh(out):
    testo = scarica(f"{GH}/doc/ghcnh-station-list.csv")
    if not testo:
        log("elenco GHCNh non raggiungibile"); return
    righe = testo.decode("latin-1").splitlines()
    log(f"GHCNh: {len(righe)} righe nell'elenco; prime: {righe[:2]}")
    scelte = []
    for r in csv.reader(righe):
        if len(r) < 6:
            continue
        try:
            lat, lon = float(r[1]), float(r[2])
        except ValueError:
            continue
        a = area_di(lat, lon)
        if a:
            scelte.append(dict(codice=r[0].strip(), nome=r[5].strip(), lat=lat, lon=lon,
                               quota=float(r[3]) if r[3].strip() else float("nan"), area=a, altro=[x.strip() for x in r[6:]]))
    log(f"GHCNh: {len(scelte)} stazioni nelle aree")
    anni = range(int(DA[:4]), time.gmtime().tm_year + 1)
    t0 = time.mktime(time.strptime(DA, "%Y-%m-%d")) - time.timezone
    intestazione_vista = False
    for s in scelte:
        T, D, V, G, K = [], [], [], [], []
        scartate = {}
        for anno in anni:
            b = None
            for url in (f"{GH}/access/by-year/{anno}/psv/GHCNh_{s['codice']}_{anno}.psv",
                        f"{GH}/access/by-year/{anno}/GHCNh_{s['codice']}_{anno}.psv"):
                b = scarica(url, tentativi=2)
                if b:
                    break
            if not b:
                continue
            lettore = csv.DictReader(io.StringIO(b.decode("utf-8", errors="replace")), delimiter="|")
            if not intestazione_vista:
                log(f"colonne GHCNh: {lettore.fieldnames}"); intestazione_vista = True
            for r in lettore:
                ff = (r.get("wind_speed") or "").strip()
                if not ff:
                    continue
                q = (r.get("wind_speed_Quality_Code") or "").strip()
                if q:
                    scartate[q] = scartate.get(q, 0) + 1; continue
                try:
                    t = time.mktime((int(r["Year"]), int(r["Month"]), int(r["Day"]), int(r["Hour"]), int(r["Minute"] or 0),
                                     0, 0, 0, 0)) - time.timezone
                except (ValueError, KeyError):
                    continue
                if t < t0:
                    continue
                dd = (r.get("wind_direction") or "").strip()
                qd = (r.get("wind_direction_Quality_Code") or "").strip()
                try:
                    d = float(dd) if dd and not qd else np.nan
                except ValueError:
                    d = np.nan
                if np.isfinite(d) and (d < 0 or d > 360):
                    d = np.nan
                gg = (r.get("wind_gust") or "").strip()
                try:
                    g = float(gg) if gg and not (r.get("wind_gust_Quality_Code") or "").strip() else np.nan
                except ValueError:
                    g = np.nan
                try:
                    v = float(ff)
                except ValueError:
                    continue
                T.append(t); D.append(d); V.append(v); G.append(g); K.append((r.get("wind_speed_Report_Type") or "").strip())
        if not T:
            log(f"  {s['codice']} {s['nome']}: nessun dato"); s["misure"] = 0; continue
        s["misure"] = salva(out, "gh", s["codice"], np.array(T, float), np.array(D, np.float32),
                            np.array(V, np.float32), np.array(G, np.float32), np.array(K))
        log(f"  {s['codice']} {s['nome']} ({s['area']}): {s['misure']} misure; scartate {scartate}")
    json.dump(scelte, open(os.path.join(out, "stazioni_gh.json"), "w"), indent=1, ensure_ascii=False)


# ---------------------------------------------------------------- ISPRA, Rete Mareografica Nazionale (porti)

SPARQL = "https://dati.isprambiente.it/sparql"
Q_STAZIONI = """PREFIX ispra-top: <https://w3id.org/italia/env/onto/top/>
PREFIX ispra-emf: <https://w3id.org/italia/env/onto/inspire-mf/>
PREFIX geo: <http://www.w3.org/2003/01/geo/wgs84_pos#>
SELECT DISTINCT str(?network) AS ?network str(?statname) AS ?statname ?stat ?lat ?long WHERE {
?stat a ispra-emf:Platform; ispra-top:name ?statname; geo:lat ?lat; geo:long ?long; ispra-emf:isPlatformOf / rdfs:label ?network .
FILTER (langMatches(lang(?network), "en")) } ORDER BY ?network ?statname"""
Q_VENTO = """PREFIX ispra-top: <https://w3id.org/italia/env/onto/top/>
PREFIX ispra-emf: <https://w3id.org/italia/env/onto/inspire-mf/>
SELECT DISTINCT str(?station) AS ?station ?stat str(?placelabel) AS ?place ?starttime ?endtime ?csvUrl WHERE {
GRAPH <https://w3id.org/italia/env/ld/rmn/> {
?obseries ispra-emf:hasObservationParameter ?parameter; ispra-top:atTime ?timeinterval; ispra-top:isPartOf ?obscolltop;
 ispra-emf:isObservationMadeBySensor ?instrument; ispra-top:hasMedia ?media .
?timeinterval ispra-top:startTime ?starttime; ispra-top:endTime ?endtime .
?parameter rdfs:label ?labelpar . ?instrument ispra-emf:isHostedBy ?stat . ?stat rdfs:label ?station .
?media ispra-top:hasDownloadURL ?csvUrl . ?obscolltop ispra-top:isPartOf ?coll .
?coll ispra-top:isPartOf <https://w3id.org/italia/env/ld/rmn/dataset>; ispra-top:isCollectionOf ?place . }
?place rdfs:label ?placelabel.
FILTER (str(?labelpar) = "Wind")
FILTER (langMatches(lang(?placelabel),"it")) FILTER (langMatches(lang(?station),"it"))
FILTER (?endtime >= "%s"^^xsd:dateTime) } ORDER BY ?csvUrl"""


def sparql(q):
    url = SPARQL + "?" + urllib.parse.urlencode({"default-graph-uri": "", "query": q, "format": "text/csv"})
    b = scarica(url, tentativi=3)
    return list(csv.DictReader(io.StringIO(b.decode("utf-8", errors="replace")))) if b else None


def rmn(out):
    st = sparql(Q_STAZIONI)
    if st is None:
        log("SPARQL di ISPRA non raggiungibile"); return
    open(os.path.join(out, "rmn_stazioni_tutte.json"), "w").write(json.dumps(st, indent=1, ensure_ascii=False))
    pos = {}
    for r in st:
        try:
            lat, lon = float(r["lat"]), float(r["long"])
        except (ValueError, KeyError):
            continue
        pos[r["stat"]] = (r["statname"], r["network"], lat, lon, area_di(lat, lon))
    log(f"RMN/RON/RMLV: {len(pos)} piattaforme; nelle aree: "
        f"{sorted(set((v[0], v[4]) for v in pos.values() if v[4]))}")
    righe = sparql(Q_VENTO % (DA + "T00:00:00Z"))
    if righe is None:
        log("SPARQL del vento non riuscito"); return
    log(f"file del vento dal {DA}: {len(righe)}; esempio {righe[:2]}")
    scelte, serie = {}, {}
    for r in righe:
        p = pos.get(r["stat"])
        if not p or not p[4]:
            continue
        codice = re.sub(r"\W+", "_", p[0].lower()).strip("_")
        scelte[codice] = dict(codice=codice, nome=p[0], lat=p[2], lon=p[3], quota=float("nan"), area=p[4], rete=p[1])
        b = scarica(r["csvUrl"], tentativi=3)
        if not b:
            continue
        testo = b.decode("utf-8", errors="replace")
        if codice not in serie:
            log(f"  {codice}: inizio del file {testo[:300]!r}")
        sep = ";" if testo.count(";") > testo.count(",") else ","
        for x in csv.DictReader(io.StringIO(testo), delimiter=sep):
            k = {kk.strip().upper(): (vv or "").strip() for kk, vv in x.items() if kk}
            ws, wd = k.get("WS", ""), k.get("WD", "")
            data = k.get("DATE") or k.get("DATA") or k.get("TIME") or k.get("DATETIME") or ""
            if not ws or not data:
                continue
            ora = k.get("TIME", "") if "DATE" in k or "DATA" in k else ""
            stringa = (data + " " + ora).strip().replace("T", " ").replace("Z", "")
            t = None
            for fm in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%Y%m%d %H%M"):
                try:
                    t = time.mktime(time.strptime(stringa[:19], fm)) - time.timezone; break
                except ValueError:
                    pass
            if t is None:
                continue
            try:
                serie.setdefault(codice, []).append((t, float(wd.replace(",", ".")) if wd else np.nan, float(ws.replace(",", "."))))
            except ValueError:
                continue
    for codice, s in serie.items():
        a = np.array(s, float)
        scelte[codice]["misure"] = salva(out, "rmn", codice, a[:, 0], a[:, 1].astype(np.float32), a[:, 2].astype(np.float32),
                                         np.full(len(a), np.nan, np.float32), np.array(["RMN"] * len(a)))
        log(f"  {codice} ({scelte[codice]['area']}): {scelte[codice]['misure']} misure")
    json.dump(list(scelte.values()), open(os.path.join(out, "stazioni_rmn.json"), "w"), indent=1, ensure_ascii=False)


# ---------------------------------------------------------------- Météo-France (orarie)

def francia(out, dipartimenti=("20", "06")):
    api = scarica("https://www.data.gouv.fr/api/1/datasets/donnees-climatologiques-de-base-horaires/")
    if not api:
        log("data.gouv non raggiungibile"); return
    j = json.loads(api)
    ris = list(j.get("resources", []))
    log(f"data.gouv: insieme {j.get('id')}, {len(ris)} risorse nella prima risposta; esempi {[r['title'] for r in ris[:5]]}")
    pagina = 1
    while True:                                   # l'API 2 elenca tutte le risorse, a pagine
        b = scarica(f"https://www.data.gouv.fr/api/2/datasets/{j.get('id')}/resources/?page={pagina}&page_size=200")
        if not b:
            break
        d = json.loads(b); nuove = d.get("data", [])
        ris += [r for r in nuove if r["id"] not in {x["id"] for x in ris}]
        if not d.get("next_page") or not nuove:
            break
        pagina += 1
    log(f"data.gouv: {len(ris)} risorse in tutto")
    for r in ris:
        if any(f"_{d}_" in r["title"] for d in dipartimenti):
            log(f"risorsa: {r['title']}  {r['url']}")
    scelte = [r for r in ris if any(f"H_{d}_" in r["title"] for d in dipartimenti)
              and r["title"].endswith(".csv.gz") and re.search(r"(latest|2024|2025|2026)", r["title"])]
    t0 = DA.replace("-", "") + "00"
    serie, info = {}, {}
    for r in scelte:
        b = scarica(r["url"])
        if not b:
            continue
        lettore = csv.DictReader(io.StringIO(gzip.decompress(b).decode("utf-8", errors="replace")), delimiter=";")
        n = 0
        for x in lettore:
            ora = x.get("AAAAMMJJHH", "")
            if ora < t0:
                continue
            try:
                lat, lon = float(x["LAT"]), float(x["LON"])
            except (ValueError, KeyError):
                continue
            a = area_di(lat, lon)
            if not a:
                continue
            ff = x.get("FF", "")
            # codici di qualità di Météo-France: 0 protetta, 1 validata, 9 filtrata (non ancora validata), 2 dubbia
            if not ff or x.get("QFF", "1") == "2":
                continue
            codice = x["NUM_POSTE"]
            info.setdefault(codice, dict(codice=codice, nome=x.get("NOM_USUEL", "").strip(), lat=lat, lon=lon,
                                         quota=float(x.get("ALTI") or "nan"), area=a))
            t = time.mktime(time.strptime(ora, "%Y%m%d%H")) - time.timezone
            dd = x.get("DD", "")
            fxi = x.get("FXI", "")
            serie.setdefault(codice, []).append((t, float(dd) if dd else np.nan, float(ff), float(fxi) if fxi else np.nan))
            n += 1
        log(f"  {r['title']}: {n} righe nelle aree")
    for codice, s in serie.items():
        a = np.array(s, float)
        info[codice]["misure"] = salva(out, "mf", codice, a[:, 0], a[:, 1].astype(np.float32), a[:, 2].astype(np.float32),
                                       a[:, 3].astype(np.float32), np.array(["MF"] * len(a)))
        log(f"  {codice} {info[codice]['nome']} ({info[codice]['area']}): {info[codice]['misure']} ore")
    json.dump(list(info.values()), open(os.path.join(out, "stazioni_mf.json"), "w"), indent=1, ensure_ascii=False)


# ---------------------------------------------------------------- sonda delle altre fonti

SONDE = {
    "meteonetwork_docs.html": "http://api.meteonetwork.it/docs/",
    "meteonetwork_stazioni_senza_chiave.txt": "https://api.meteonetwork.it/v3/stations",
    "ispra_lod.html": "https://dati.isprambiente.it/?p=33",
    "mareografico.html": "https://www.mareografico.it/",
    "sir_toscana.html": "https://www.sir.toscana.it/consistenza-rete",
}


def sonda(out):
    for nome, url in SONDE.items():
        b = scarica(url, tentativi=2)
        if b is None:
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "lens-banco/1.0"})
                urllib.request.urlopen(req, timeout=60)
            except urllib.error.HTTPError as e:
                b = f"HTTP {e.code}\n".encode() + (e.read() or b"")
            except Exception as e:
                b = f"{type(e).__name__} {e}".encode()
        open(os.path.join(out, nome), "wb").write(b or b"")
        log(f"{nome}: {len(b or b'')} byte")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--isd", action="store_true")
    ap.add_argument("--francia", action="store_true")
    ap.add_argument("--ghcnh", action="store_true")
    ap.add_argument("--rmn", action="store_true")
    ap.add_argument("--sonda", action="store_true")
    ap.add_argument("--out", default="stazioni")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    try:
        if a.sonda:
            sonda(a.out)
        if a.isd:
            isd(a.out)
        if a.ghcnh:
            ghcnh(a.out)
        if a.rmn:
            rmn(a.out)
        if a.francia:
            francia(a.out)
    finally:
        open(os.path.join(a.out, f"registro_{int(time.time())}.txt"), "w").write("\n".join(LOG))
