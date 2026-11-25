"""Stazioni a terra con il vento vicino alle aree del banco: aeroporti e stazioni sinottiche (NOAA ISD),
Météo-France (Corsica e Alpi Marittime), più una sonda per le altre fonti.

  python stazioni.py --isd --out stazioni          aeroporti e stazioni sinottiche (anche fari e capi dell'Aeronautica)
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


# ---------------------------------------------------------------- Météo-France (orarie)

def francia(out, dipartimenti=("20", "06")):
    api = scarica("https://www.data.gouv.fr/api/1/datasets/donnees-climatologiques-de-base-horaires/")
    if not api:
        log("data.gouv non raggiungibile"); return
    ris = json.loads(api)["resources"]
    for r in ris:
        if any(f"H_{d}_" in r["title"] for d in dipartimenti):
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
    ap.add_argument("--sonda", action="store_true")
    ap.add_argument("--out", default="stazioni")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    try:
        if a.sonda:
            sonda(a.out)
        if a.isd:
            isd(a.out)
        if a.francia:
            francia(a.out)
    finally:
        open(os.path.join(a.out, f"registro_{int(time.time())}.txt"), "w").write("\n".join(LOG))
