"""Boe e stazioni in mare con il vento, dal prodotto in situ di Copernicus Marine (Mediterraneo), sulle aree del banco.

  python insitu.py --indici --out insitu                    scarica gli indici e dice quali piattaforme cadono nelle aree
  python insitu.py --scarica --out insitu [--da 2024-06-01] scarica le serie e tiene vento e raffica

Prodotto INSITU_MED_PHYBGCWAV_DISCRETE_MYNRT_013_035 (dataset cmems_obs-ins_med_phybgcwav_mynrt_na_irr): un file NetCDF
per piattaforma (cartella history), con WSPD (vento), WDIR, GSPD (raffica) e le loro bandiere di qualità.
Uscita: insitu_<piattaforma>.npz con t (s UTC), lat, lon, vento (m/s), direzione, raffica, quota del sensore se c'è.
Credenziali: CMEMS_USER e CMEMS_PASS (segreti del deposito).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from costa import AREE

DS = "cmems_obs-ins_med_phybgcwav_mynrt_na_irr"
LOG = []


def log(x):
    LOG.append(str(x)); print(x, flush=True)


def login():
    import copernicusmarine as cm
    u, p = os.environ.get("CMEMS_USER"), os.environ.get("CMEMS_PASS")
    return dict(username=u, password=p) if u and p else {}


def indici(out):
    import copernicusmarine as cm
    cred = login()
    for nome in ("index_platform.txt", "index_history.txt"):
        try:
            cm.get(dataset_id=DS, filter=f"*{nome}", output_directory=out, no_directories=True, **cred)
        except Exception as e:
            log(f"{nome}: {type(e).__name__} {e}")
    trovati = glob.glob(os.path.join(out, "**", "index_history.txt"), recursive=True)
    if not trovati:
        log("index_history.txt non trovato"); return []
    righe = [r for r in open(trovati[0], encoding="utf-8", errors="replace") if not r.startswith("#")]
    log(f"index_history: {len(righe)} righe; esempio: {righe[0][:300] if righe else ''}")
    scelti = []
    for r in righe:
        c = [x.strip() for x in r.split(",")]
        # catalog_id,file_name,geospatial_lat_min,geospatial_lat_max,geospatial_lon_min,geospatial_lon_max,
        # time_coverage_start,time_coverage_end,provider,date_update,data_mode,parameters
        try:
            f = c[1]; la0, la1, lo0, lo1 = map(float, c[2:6]); fine = c[7]; par = c[-1]
        except Exception:
            continue
        if "WSPD" not in par or fine < "2024-06":
            continue
        for area, (s, n, w, e) in AREE.items():
            if la1 >= s - 0.3 and la0 <= n + 0.3 and lo1 >= w - 0.3 and lo0 <= e + 0.3 and (la1 - la0) < 0.5:
                scelti.append({"file": f, "area": area, "lat": [la0, la1], "lon": [lo0, lo1], "fine": fine, "parametri": par})
    log(f"{len(scelti)} serie con il vento nelle aree")
    for x in scelti:
        log(json.dumps(x))
    json.dump(scelti, open(os.path.join(out, "scelti.json"), "w"), indent=1)
    return scelti


def scarica(out, da):
    import copernicusmarine as cm
    import xarray as xr
    scelti = json.load(open(os.path.join(out, "scelti.json")))
    cred = login()
    t0 = np.datetime64(da)
    for x in scelti:
        nome = os.path.basename(x["file"])
        try:
            cm.get(dataset_id=DS, filter=f"*{nome}", output_directory=os.path.join(out, "nc"), no_directories=True, **cred)
            p = glob.glob(os.path.join(out, "nc", "**", nome), recursive=True)[0]
            d = xr.open_dataset(p, decode_times=True)
            t = d["TIME"].values
            m = t >= t0

            def var(v, qc=True):
                if v not in d:
                    return None
                a = d[v].values
                if a.ndim == 2:          # (tempo, profondità/quota): il primo livello valido
                    a = np.array([r[np.isfinite(r)][0] if np.isfinite(r).any() else np.nan for r in a])
                q = d.get(v + "_QC")
                if qc and q is not None:
                    qv = q.values
                    if qv.ndim == 2:
                        qv = qv[:, 0]
                    a = np.where(np.isin(qv, [1, 2]), a, np.nan)
                return a[m].astype("f4")
            ws = var("WSPD")
            if ws is None or np.isfinite(ws).sum() < 100:
                log(f"{nome}: vento assente o troppo poco"); continue
            lat = d["LATITUDE"].values; lon = d["LONGITUDE"].values
            lat = np.broadcast_to(lat, t.shape)[m] if lat.size > 1 else np.full(m.sum(), float(lat))
            lon = np.broadcast_to(lon, t.shape)[m] if lon.size > 1 else np.full(m.sum(), float(lon))
            altezza = None
            for k in ("HEIGHT", "DEPH"):
                if k in d:
                    altezza = float(np.nanmedian(d[k].values)); break
            np.savez_compressed(os.path.join(out, f"insitu_{nome.replace('.nc', '')}.npz"),
                                t=(t[m].astype("datetime64[s]").astype("i8")), lat=lat.astype("f4"), lon=lon.astype("f4"),
                                vento=ws, direzione=var("WDIR") if var("WDIR") is not None else np.full(m.sum(), np.nan, "f4"),
                                raffica=var("GSPD") if var("GSPD") is not None else np.full(m.sum(), np.nan, "f4"),
                                meta=np.array(json.dumps({**x, "altezza": altezza, "attributi": {k: str(v) for k, v in d.attrs.items()
                                                                                               if k in ("platform_name", "platform_code", "institution", "wmo_platform_code")}})))
            log(f"{nome} ({x['area']}): {np.isfinite(ws).sum()} misure dal {da}, quota {altezza}")
            d.close()
        except Exception as e:
            log(f"{nome}: {type(e).__name__} {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indici", action="store_true"); ap.add_argument("--scarica", action="store_true")
    ap.add_argument("--out", default="insitu"); ap.add_argument("--da", default="2024-06-01")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.indici or a.scarica:
        indici(a.out)
    if a.scarica:
        scarica(a.out, a.da)
    open(os.path.join(a.out, "log_insitu.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
