"""Le corse ECMWF dei dati aperti per i passaggi del radar: archivio su cloud pubblico, dal 2024.

Il radar dei Sentinel-1 passa verso le 05 e le 17 UTC. La corsa più recente prima del passaggio è quella delle 00 o
delle 12, con 4-6 ore di anticipo: bastano i passi a 3 e 6 ore. Per ogni giorno del mese si scaricano, dal file
globale a 0,25 gradi, solo il vento a 10 m (10u, 10v) con una richiesta a pezzi guidata dal file .index; si ritaglia
sul riquadro delle aree di Lens e si scrive ecmwf_<AAAAMMGGHH>.grb2, come i file dell'archivio di shift-grib.

  python ecmwf_aperto.py --mese 2025-07 --out ecmwf [--box 35,45.5,7.5,16.5] [--passi 3,6]

Dati: ECMWF open data, CC BY 4.0 (copia pubblica su AWS, registro Open Data).
Il ritaglio è lo stesso di shift-grib (tools/build_gribs.py, crop_message).
"""
import argparse
import calendar
import datetime as dt
import json
import os
import sys
import time
import urllib.request

import eccodes as ec
import numpy as np

FONTI = [
    "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com",
    "https://ai4edataeuwest.blob.core.windows.net/ecmwf",
    "https://data.ecmwf.int/forecasts",
]


def get(url, rng=None, tries=3):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
            time.sleep(2 + 3 * k)
        except Exception:
            time.sleep(2 + 3 * k)
    return None


def nome(root, t, step, ext):
    return f"{root}/{t:%Y%m%d}/{t:%H}z/ifs/0p25/oper/{t:%Y%m%d%H}0000-{step}h-oper-fc.{ext}"


def crop_message(h, box):
    if ec.codes_get(h, "gridType") != "regular_ll":
        return None
    ni, nj = ec.codes_get(h, "Ni"), ec.codes_get(h, "Nj")
    la1 = ec.codes_get(h, "latitudeOfFirstGridPointInDegrees"); la2 = ec.codes_get(h, "latitudeOfLastGridPointInDegrees")
    lo1 = ec.codes_get(h, "longitudeOfFirstGridPointInDegrees")
    di = ec.codes_get(h, "iDirectionIncrementInDegrees"); dj = ec.codes_get(h, "jDirectionIncrementInDegrees")
    if ec.codes_get(h, "iScansNegatively") or ec.codes_get(h, "jPointsAreConsecutive"):
        return None
    lats = la1 + np.arange(nj) * (dj if la2 > la1 else -dj)
    lons = ((lo1 + np.arange(ni) * di + 180.0) % 360.0) - 180.0
    latmin, latmax, lonmin, lonmax = box
    jj = np.where((lats >= latmin - 1e-6) & (lats <= latmax + 1e-6))[0]
    order = np.argsort(lons, kind="stable")
    ii = order[(lons[order] >= lonmin - 1e-6) & (lons[order] <= lonmax + 1e-6)]
    if len(jj) < 2 or len(ii) < 2:
        return None
    vals = ec.codes_get_values(h).reshape(nj, ni)
    sub = vals[np.ix_(jj, ii)]
    c = ec.codes_clone(h)
    try:
        ec.codes_set(c, "packingType", "grid_simple")
        ec.codes_set(c, "Ni", len(ii)); ec.codes_set(c, "Nj", len(jj))
        ec.codes_set(c, "latitudeOfFirstGridPointInDegrees", float(lats[jj[0]]))
        ec.codes_set(c, "latitudeOfLastGridPointInDegrees", float(lats[jj[-1]]))
        ec.codes_set(c, "longitudeOfFirstGridPointInDegrees", float(lons[ii[0]] % 360.0))
        ec.codes_set(c, "longitudeOfLastGridPointInDegrees", float(lons[ii[-1]] % 360.0))
        ec.codes_set(c, "bitsPerValue", 12)
        ec.codes_set_values(c, sub.flatten())
        return ec.codes_get_message(c)
    finally:
        ec.codes_release(c)


def crop_bytes(data, box):
    out, pos = [], 0
    while True:
        i = data.find(b"GRIB", pos)
        if i < 0:
            return out
        h = ec.codes_new_from_message(data[i:])
        try:
            n = ec.codes_get(h, "totalLength")
            m = crop_message(h, box)
            if m:
                out.append(m)
        finally:
            ec.codes_release(h)
        pos = i + max(int(n), 4)


def corsa(t, passi, box, log):
    for root in FONTI:
        msgs = []
        for step in passi:
            idx = get(nome(root, t, step, "index"))
            if not idx:
                msgs = None
                break
            for line in idx.decode().splitlines():
                rec = json.loads(line)
                if rec.get("param") in ("10u", "10v") and rec.get("levtype", "sfc") == "sfc":
                    a, n = int(rec["_offset"]), int(rec["_length"])
                    raw = get(nome(root, t, step, "grib2"), (a, a + n - 1))
                    if raw is None:
                        msgs = None
                        break
                    msgs += crop_bytes(raw, box)
            if msgs is None:
                break
        if msgs:
            return msgs, root
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mese", required=True)
    ap.add_argument("--out", default="ecmwf")
    ap.add_argument("--box", default="35,45.5,7.5,16.5")
    ap.add_argument("--passi", default="3,6")
    a = ap.parse_args()
    box = [float(x) for x in a.box.split(",")]
    passi = [int(x) for x in a.passi.split(",")]
    os.makedirs(a.out, exist_ok=True)
    y, m = [int(x) for x in a.mese.split("-")]
    log = []
    fatte = 0
    for d in range(1, calendar.monthrange(y, m)[1] + 1):
        for h in (0, 12):
            t = dt.datetime(y, m, d, h)
            dst = os.path.join(a.out, f"ecmwf_{t:%Y%m%d%H}.grb2")
            if os.path.exists(dst):
                continue
            msgs, root = corsa(t, passi, box, log)
            if not msgs:
                log.append(f"{t:%Y%m%d%H}: non trovata")
                continue
            open(dst, "wb").write(b"".join(msgs))
            fatte += 1
            log.append(f"{t:%Y%m%d%H}: {len(msgs)} campi da {root.split('/')[2]}, {os.path.getsize(dst)} byte")
            print(log[-1], flush=True)
    open(os.path.join(a.out, f"log_ecmwf_{a.mese}.txt"), "w").write("\n".join(log) + "\n")
    print(f"{a.mese}: {fatte} corse")
    if fatte == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
