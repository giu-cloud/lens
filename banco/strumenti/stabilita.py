"""Stabilità e contrasto termico della scena, dai campi aperti di ECMWF (2 m, mare, 925 e 850 hPa, 100 m, raffica).

  python3 stabilita.py <cartella stab_*.grb2> <cartella del metro (campi_*.npz)> <rilievo .rgr> <uscita> [cartella ecmwf_*.grb2]
(il vento a 10 m per taglio e raffica viene dai file del vento, ecmwf_*.grb2, della stessa corsa)

Per ogni area e passaggio, con la stessa corsa e lo stesso anticipo del vento (la corsa delle 00 o 12 prima del
passaggio, interpolata tra i passi a 3 e 6 ore), sui punti della griglia a 0,25° nel riquadro dell'area (± 0,25°):
  sst       temperatura della superficie del mare (skt sui punti di mare), °C
  dTam      2 m sul mare − mare (K): negativo = aria più fredda del mare, instabile; positivo = stabile
  dT925     θ a 925 hPa − θ del mare (K): quanto è forte lo strato stabile sopra il mare (inversione bassa)
  dT850     θ a 850 hPa − θ del mare (K)
  terra_mare  2 m sulla terra (entro ± 0,5°) − mare (K): motore delle brezze
  taglio    vento a 100 m ÷ vento a 10 m sul mare (oltre 1,25 di solito aria stabile)
  raff      raffica ÷ vento a 10 m sul mare (alta con aria instabile)
  v925, d925  forza (m/s) e direzione di provenienza (gradi) a 925 hPa; v850, d850 lo stesso a 850
  umid      2 m: differenza tra temperatura e rugiada sul mare (K)
"""
import datetime as dt
import os
import sys

import eccodes as ec
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rilievo import Rilievo

AREE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]
NOMI = ["sst", "dTam", "dT925", "dT850", "terra_mare", "taglio", "raff", "v925", "d925", "v850", "d850", "umid"]


def leggi(p):
    campi = {}
    lat = lon = None
    with open(p, "rb") as f:
        while True:
            h = ec.codes_grib_new_from_file(f)
            if h is None:
                break
            try:
                sn = ec.codes_get(h, "shortName"); lv = ec.codes_get(h, "level"); st = ec.codes_get(h, "step")
                tl = ec.codes_get(h, "typeOfLevel")
                ni, nj = ec.codes_get(h, "Ni"), ec.codes_get(h, "Nj")
                v = ec.codes_get_values(h).reshape(nj, ni)
                if lat is None:
                    la1 = ec.codes_get(h, "latitudeOfFirstGridPointInDegrees"); la2 = ec.codes_get(h, "latitudeOfLastGridPointInDegrees")
                    lo1 = ec.codes_get(h, "longitudeOfFirstGridPointInDegrees"); di = ec.codes_get(h, "iDirectionIncrementInDegrees")
                    lat = np.linspace(la1, la2, nj); lon = lo1 + np.arange(ni) * di
                    lon = np.where(lon > 180, lon - 360, lon)
                chiave = sn if tl in ("surface", "meanSea") else (f"{sn}@{lv}" if tl == "isobaricInhPa" else f"{sn}{lv}" if sn in ("u", "v") else sn)
                chiave = {"100u": "u100", "100v": "v100", "max_i10fg": "10fg", "i10fg": "10fg"}.get(chiave, chiave)
                campi[(chiave, st)] = v
            finally:
                ec.codes_release(h)
    return lat, lon, campi


def corsa_per(t, disp):
    d = dt.datetime.fromtimestamp(int(t), dt.timezone.utc).replace(tzinfo=None)
    for back in range(0, 4):
        base = d.replace(minute=0, second=0, microsecond=0)
        r = base.replace(hour=(base.hour // 12) * 12) - dt.timedelta(hours=12 * back)
        lead = (d - r).total_seconds() / 3600
        if lead > 18:
            break
        n = f"stab_{r:%Y%m%d%H}.grb2"
        if n in disp:
            return n, lead
    return None, None


def theta(T, p):
    return T * (1000.0 / p) ** 0.2857


def main():
    cart, metro, ril_p, usc = sys.argv[1:5]
    vento = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(usc, exist_ok=True)
    ril = Rilievo(ril_p)
    disp = set(os.listdir(cart))
    cache = {}
    for area in AREE:
        z = np.load(os.path.join(metro, f"campi_{area}.npz"))
        t = z["t"]
        la0, la1, lo0, lo1 = z["lat"].min(), z["lat"].max(), z["lon"].min(), z["lon"].max()
        out = {k: np.full(len(t), np.nan, dtype="f4") for k in NOMI}
        maschere = None
        for s, ts in enumerate(t):
            n, lead = corsa_per(ts, disp)
            if n is None:
                continue
            if n not in cache:
                if len(cache) > 40:
                    cache.clear()
                lat_, lon_, C_ = leggi(os.path.join(cart, n))
                pv = os.path.join(vento, n.replace("stab_", "ecmwf_")) if vento else None
                if pv and os.path.exists(pv):
                    C_.update(leggi(pv)[2])
                cache[n] = (lat_, lon_, C_)
            lat, lon, C = cache[n]
            if maschere is None:
                LA, LO = np.meshgrid(lat, lon, indexing="ij")
                terra = np.zeros(LA.shape, bool)
                for j in range(LA.shape[0]):
                    for i in range(LA.shape[1]):
                        if la0 - 0.6 <= LA[j, i] <= la1 + 0.6 and lo0 - 0.6 <= LO[j, i] <= lo1 + 0.6:
                            _, _, hh = ril.griglia(0, LA[j, i] - 0.04, LA[j, i] + 0.04, LO[j, i] - 0.05, LO[j, i] + 0.05)
                            terra[j, i] = (hh > 0).mean() > 0.5
                vicino = (LA >= la0 - 0.25) & (LA <= la1 + 0.25) & (LO >= lo0 - 0.25) & (LO <= lo1 + 0.25)
                largo = (LA >= la0 - 0.5) & (LA <= la1 + 0.5) & (LO >= lo0 - 0.5) & (LO <= lo1 + 0.5)
                mare = vicino & ~terra
                if mare.sum() < 2:
                    mare = largo & ~terra
                maschere = (mare, largo & terra)
            mare, terra_v = maschere
            w = (lead - 3.0) / 3.0
            w = min(max(w, 0.0), 1.0) if 3 <= lead <= 6 else None

            def campo(k):
                a, b = C.get((k, 3)), C.get((k, 6))
                if a is None or b is None:
                    return None
                if w is None:
                    return a if abs(lead - 3) < abs(lead - 6) else b
                return (1 - w) * a + w * b

            def media(k, m):
                f = campo(k)
                return float(np.mean(f[m])) if f is not None and m.sum() else np.nan

            sst = media("skt", mare); t2 = media("2t", mare); d2 = media("2d", mare); msl = media("msl", mare) / 100.0
            ths = theta(sst, msl) if np.isfinite(msl) else np.nan
            out["sst"][s] = sst - 273.15
            out["dTam"][s] = t2 - sst
            out["dT925"][s] = theta(media("t@925", mare), 925.0) - ths
            out["dT850"][s] = theta(media("t@850", mare), 850.0) - ths
            out["terra_mare"][s] = media("2t", terra_v) - sst if terra_v.sum() else np.nan
            out["umid"][s] = t2 - d2
            u10 = campo("10u"); v10 = campo("10v")
            u100 = campo("u100"); v100 = campo("v100"); fg = campo("10fg")
            if u100 is not None and v100 is not None and u10 is not None:
                s10 = np.hypot(u10, v10)[mare]; s100 = np.hypot(u100, v100)[mare]
                out["taglio"][s] = float(np.mean(s100) / max(np.mean(s10), 0.5))
            if fg is not None and u10 is not None:
                out["raff"][s] = float(np.mean(fg[mare]) / max(np.mean(np.hypot(u10, v10)[mare]), 0.5))
            for lv in ("925", "850"):
                uu = campo(f"u@{lv}"); vv = campo(f"v@{lv}")
                if uu is not None and vv is not None:
                    um, vm = float(np.mean(uu[mare | terra_v])), float(np.mean(vv[mare | terra_v]))
                    out[f"v{lv}"][s] = np.hypot(um, vm)
                    out[f"d{lv}"][s] = (np.degrees(np.arctan2(-um, -vm)) + 360) % 360
        np.savez(os.path.join(usc, f"stab_{area}.npz"), t=t, **out)
        ok = np.isfinite(out["dT925"])
        print(f"{area}: {ok.sum()} passaggi su {len(t)}; dT925 mediana {np.nanmedian(out['dT925']):.1f} K, "
              f"dTam {np.nanmedian(out['dTam']):+.1f} K", flush=True)


if __name__ == "__main__":
    main()
