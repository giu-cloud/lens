"""Storia di ICON-2I (ItaliaMeteo-Arpae) dall'archivio aperto di Open-Meteo su AWS, ritagliata sulle aree del banco.

  python icon2i_storia.py --aree palermo,egadi --da 2025-04-10 --a 2026-11-25 --out icon2i [--ore 4,5,6,16,17,18]

Open-Meteo tiene ICON-2I in file .om a blocchi di 96 ore (chunk_<k>.om, ore da k*96 dall'epoca), con la griglia
intera (33,7-48,9 N, 3-22 E). Si scarica un blocco alla volta, se ne tengono solo le aree e le ore che servono e si cancella: vento a 10 m
(forza e direzione) nelle ore dei passaggi del radar (di default 4-6 e 16-18 UTC).
È la serie "storica" di Open-Meteo: per ogni ora le prime ore della corsa più recente (anticipo di poche ore).
Uscita: icon2i_<area>.npz con lat, lon, t (secondi UTC), u, v (m/s, float16) e un registro.
Dati: Open-Meteo (CC BY 4.0), modello ItaliaMeteo-Arpae ICON-2I.
"""
import argparse
import datetime as dt
import json
import os
import sys
import time

import numpy as np
from omfiles import OmFileReader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from costa import AREE

BASE = "https://openmeteo.s3.amazonaws.com/data/italia_meteo_arpae_icon_2i"
BBOX = (33.7, 48.9, 3.0, 22.0)
ORA = 3600
BLOCCO = 96


def apri(fs, var, k, cartella="/tmp/om"):
    """Scarica il blocco intero (15-30 MB) e lo apre; a fine uso il file si cancella."""
    import requests
    os.makedirs(cartella, exist_ok=True)
    dst = os.path.join(cartella, f"{var}_{k}.om")
    for tent in range(4):
        try:
            r = requests.get(f"{BASE}/{var}/chunk_{k}.om", timeout=120)
            if r.status_code == 404:
                raise FileNotFoundError(f"404 {var} {k}")
            r.raise_for_status()
            open(dst, "wb").write(r.content)
            return OmFileReader.from_path(dst)
        except FileNotFoundError:
            raise
        except Exception as e:
            err = e; time.sleep(5 * (tent + 1))
    raise RuntimeError(f"{var} {k}: {err}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aree", required=True)
    ap.add_argument("--da", default="2025-04-10")
    ap.add_argument("--a", default=dt.date.today().isoformat())
    ap.add_argument("--out", default="icon2i")
    ap.add_argument("--ore", default="4,5,6,16,17,18")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    aree = [x for x in a.aree.split(",") if x in AREE]
    ore = {int(x) for x in a.ore.split(",")}
    fs = None
    t0 = int(dt.datetime.fromisoformat(a.da).replace(tzinfo=dt.timezone.utc).timestamp())
    t1 = int(dt.datetime.fromisoformat(a.a).replace(tzinfo=dt.timezone.utc).timestamp())
    k0, k1 = t0 // (BLOCCO * ORA), t1 // (BLOCCO * ORA)
    log = []
    acc = {ar: {"t": [], "u": [], "v": []} for ar in aree}
    griglia = None
    for k in range(k0, k1 + 1):
        try:
            rs = apri(fs, "wind_speed_10m", k); rd = apri(fs, "wind_direction_10m", k)
        except Exception as e:
            log.append(f"chunk {k}: non c'è ({type(e).__name__}: {e})")
            continue
        sh = rs.shape
        if griglia is None:
            ny, nx = sh[0], sh[1]
            lat = np.linspace(BBOX[0], BBOX[1], ny); lon = np.linspace(BBOX[2], BBOX[3], nx)
            griglia = (lat, lon)
            log.append(f"forma {sh}, blocchi {rs.chunks}, passo lat {lat[1]-lat[0]:.5f} lon {lon[1]-lon[0]:.5f}")
            # controllo delle unità: il punto 38,2 N 13,3 E alle 00 del 1/5/2025 (l'API dava 2,5 km/h, da 183°)
        lat, lon = griglia
        nt = sh[2]
        ts = (k * BLOCCO + np.arange(nt)) * ORA
        sel = np.array([((tt // ORA) % 24) in ore for tt in ts])
        if not sel.any():
            continue
        i0, i1 = np.nonzero(sel)[0][[0, -1]]
        for ar in aree:
            s, n, w, e = AREE[ar]
            ja = int(np.searchsorted(lat, s - 0.05)); jb = int(np.searchsorted(lat, n + 0.05))
            ia = int(np.searchsorted(lon, w - 0.05)); ib = int(np.searchsorted(lon, e + 0.05))
            for tent in range(3):
                try:
                    S = rs[ja:jb, ia:ib, i0:i1 + 1].astype("f4")
                    D = rd[ja:jb, ia:ib, i0:i1 + 1].astype("f4")
                    break
                except Exception as ex:
                    log.append(f"chunk {k} {ar}: errore {ex}, riprovo"); time.sleep(5)
            else:
                continue
            ss = sel[i0:i1 + 1]
            S = S[:, :, ss]; D = D[:, :, ss]
            u = -S * np.sin(np.radians(D)); v = -S * np.cos(np.radians(D))
            acc[ar]["t"].append(ts[i0:i1 + 1][ss]); acc[ar]["u"].append(np.moveaxis(u, 2, 0)); acc[ar]["v"].append(np.moveaxis(v, 2, 0))
            if "lat" not in acc[ar]:
                acc[ar]["lat"] = lat[ja:jb]; acc[ar]["lon"] = lon[ia:ib]
        if (k - k0) % 20 == 0:
            log.append(f"chunk {k} ({dt.datetime.fromtimestamp(k * BLOCCO * ORA, dt.timezone.utc):%Y-%m-%d}) fatto")
            print(log[-1], flush=True)
        # controllo delle unità sul punto di prova
        tp = int(dt.datetime(2025, 5, 1, tzinfo=dt.timezone.utc).timestamp())
        if ts[0] <= tp <= ts[-1]:
            jp = int(np.argmin(abs(lat - 38.2))); ip = int(np.argmin(abs(lon - 13.3))); ip_t = int((tp - ts[0]) // ORA)
            log.append(f"punto di prova 38,2 13,3 alle 00 del 1/5/2025: forza {float(rs[jp, ip, ip_t]):.2f}, "
                       f"direzione {float(rd[jp, ip, ip_t]):.0f} (l'API dava 2,5 km/h da 183°)")
        rs.close(); rd.close()
        for f in os.listdir("/tmp/om"):
            os.remove(os.path.join("/tmp/om", f))
    for ar in aree:
        d = acc[ar]
        if not d["t"]:
            log.append(f"{ar}: niente"); continue
        np.savez_compressed(os.path.join(a.out, f"icon2i_{ar}.npz"), lat=d["lat"], lon=d["lon"],
                            t=np.concatenate(d["t"]), u=np.concatenate(d["u"]).astype("f2"),
                            v=np.concatenate(d["v"]).astype("f2"),
                            meta=np.array(json.dumps({"fonte": "Open-Meteo, ItaliaMeteo-Arpae ICON-2I", "unita": "come nel file"})))
        log.append(f"{ar}: {sum(len(x) for x in d['t'])} ore, griglia {len(d['lat'])}x{len(d['lon'])}")
    open(os.path.join(a.out, f"log_icon2i_{a.aree.replace(',', '_')[:60]}.txt"), "w").write("\n".join(log) + "\n")
    print("\n".join(log[-20:]))


if __name__ == "__main__":
    main()
