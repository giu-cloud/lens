"""Il radar contro ECMWF come lo legge Shift (GRIB 0,25 gradi dell'archivio): lettura normale, lettura dai soli punti di
mare, e Lens sopra ciascuna. Per distanza dalla terra: rapporto mediano radar/modello, errore medio e quadratico in nodi.

  python3 costa_ecmwf.py <dati> <atlanti> <ecmwf> <uscita> area [area ...]
"""
import datetime as dt
import json
import os
import subprocess
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import costa_analisi as ca

S = os.path.dirname(os.path.abspath(__file__))
KN = 1.943844
BINS = [(0, 1.5), (1.5, 3), (3, 5), (5, 10), (10, 999)]
NOMI = ["cmems", "shift", "mare", "lens", "lens_mare"]


def run_per(t, disponibili):
    """La corsa più recente partita prima della scena, entro 18 ore."""
    d = dt.datetime.utcfromtimestamp(t)
    for back in range(0, 4):
        base = d.replace(minute=0, second=0, microsecond=0)
        h = (base.hour // 12) * 12
        r = base.replace(hour=h) - dt.timedelta(hours=12 * back)
        lead = (d - r).total_seconds() / 3600
        if lead > 18:
            break
        n = f"ecmwf_{r:%Y%m%d%H}.grb2"
        if n in disponibili:
            return n, lead
    return None, None


def main():
    dati, atl, ecm, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    disp = set(os.listdir(ecm))
    risultati = {}
    for area in sys.argv[5:]:
        lat, lon, scene, _ = ca.carica(dati, area)
        sc = []
        for s in scene:
            n, lead = run_per(s["t"], disp)
            if n:
                sc.append((s, n, lead))
        if not sc:
            print(f"{area}: nessuna scena con l'archivio"); continue
        ny, nx = len(lat), len(lon)
        ib = os.path.join(out, f"sar2_{area}.bin"); ob = os.path.join(out, f"sho_{area}.bin")
        with open(ib, "wb") as f:
            f.write(b"SAR2"); f.write(np.array([ny, nx], "<i4").tobytes())
            f.write(np.array([lat[0], lat[1] - lat[0], lon[0], lon[1] - lon[0]], "<f8").tobytes())
            f.write(np.array([len(sc)], "<i4").tobytes())
            for s, n, _ in sc:
                f.write(np.array([s["t"] * 1000], "<i8").tobytes())
                b = n.encode(); f.write(np.array([len(b)], "<i4").tobytes()); f.write(b)
        atlante = os.path.join(atl, f"{area}.lna")
        cp = f"{S}/shiftsar.jar:{os.environ.get('UT', S + '/../ut')}:/tmp/kotlinc/lib/kotlin-stdlib.jar"
        props = [f"-D{k}={v}" for k, v in json.loads(os.environ.get("LENS_PROPS", "{}")).items()]
        sea = [f"-D{k}={v}" for k, v in json.loads(os.environ.get("SEA_PROPS", "{}")).items()]
        usa = atlante if os.path.exists(atlante) and not os.environ.get("NOLENS") else "-"
        subprocess.run(["java", "-Xmx4g", *props, *sea, "-cp", cp, "ShiftSarKt", usa, ib, ob, ecm], check=True)
        raw = np.fromfile(ob, dtype="<f4", offset=16).reshape(len(sc), 8, ny, nx).astype("float64")
        SR = np.array([np.hypot(s["su"], s["sv"]) for s, _, _ in sc])
        CM = np.array([np.hypot(s["mu"], s["mv"]) for s, _, _ in sc])
        campi = {"cmems": CM, "shift": np.hypot(raw[:, 0], raw[:, 1]), "mare": np.hypot(raw[:, 2], raw[:, 3]),
                 "lens": np.hypot(raw[:, 4], raw[:, 5]), "lens_mare": np.hypot(raw[:, 6], raw[:, 7])}
        # anche la direzione: errore medio assoluto in gradi
        dirs = {"radar": np.degrees(np.arctan2(np.array([s["su"] for s, _, _ in sc]), np.array([s["sv"] for s, _, _ in sc]))),
                "shift": np.degrees(np.arctan2(raw[:, 0], raw[:, 1])), "mare": np.degrees(np.arctan2(raw[:, 2], raw[:, 3])),
                "lens": np.degrees(np.arctan2(raw[:, 4], raw[:, 5])), "lens_mare": np.degrees(np.arctan2(raw[:, 6], raw[:, 7]))}
        mai = ~np.isfinite(SR).any(axis=0)
        kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
        dist = ndimage.distance_transform_edt(~mai, sampling=(kmy, kmx))
        D = np.broadcast_to(dist, SR.shape)
        lensok = not os.environ.get("NOLENS")
        nomi = NOMI if lensok else NOMI[:3]
        base = np.isfinite(SR) & np.all([np.isfinite(campi[n]) for n in nomi], axis=0) & (campi["shift"] > 3.0)
        leads = np.array([l for _, _, l in sc])
        res = {"scene": len(sc), "lead_medio_h": float(leads.mean()), "celle": int(base.sum()), "bin": {}}
        print(f"\n== {area}: {len(sc)} scene con l'archivio (anticipo medio {leads.mean():.1f} h), {base.sum()} celle con Shift oltre 6 nodi")
        print(f"   {'km dalla terra':<14}" + "".join(f"{n:>22}" for n in nomi) + "     celle")
        for lo, hi in BINS:
            m = base & (D >= lo) & (D < hi)
            if m.sum() < 50:
                continue
            r = {}
            for n in nomi:
                x = campi[n][m]; y = SR[m]
                r[n] = dict(rapporto=float(np.median(y / x)), errore_kn=float(np.mean(x - y) * KN), rmse_kn=float(np.sqrt(np.mean((x - y) ** 2)) * KN))
                if n in dirs:
                    dd = (dirs[n][m] - dirs["radar"][m] + 540) % 360 - 180
                    r[n]["dir_gradi"] = float(np.mean(np.abs(dd)))
            r["celle"] = int(m.sum())
            r["radar_kn"] = float(np.mean(SR[m]) * KN)
            res["bin"][f"{lo}-{hi}"] = r
            print(f"   {lo:>4}-{hi:<4} km   " + "".join(f"  {r[n]['rapporto']:.2f} {r[n]['errore_kn']:+5.1f} {r[n]['rmse_kn']:4.1f}kn" for n in nomi) + f"  {m.sum():>8}")
        print("   (per colonna: radar/modello mediano, errore medio del modello in nodi, errore quadratico in nodi)")
        risultati[area] = res
        np.savez_compressed(os.path.join(out, f"campi_{area}.npz"), lat=lat, lon=lon, dist=dist, t=np.array([s["t"] for s, _, _ in sc]),
                            radar=SR.astype("f4"), **{k: v.astype("f4") for k, v in campi.items()})
    with open(os.path.join(out, f"ecmwf_shift_{'_'.join(sys.argv[5:]) if len(sys.argv) < 8 else 'tutte'}.json"), "w") as f:
        json.dump(risultati, f, indent=1)


if __name__ == "__main__":
    main()
