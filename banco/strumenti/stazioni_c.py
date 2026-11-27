"""C (Lens Fast con il livello del radar) contro le stazioni a terra tarate dal radar, a tutte le ore (27/11 sera).

  SUFFISSO=_s025 PROPS="-Dshift.lens.relief=0.25" [KM=10] [MAX=400] python3 stazioni_c.py <area> [<area> ...]

Parte dai confronti di stazioni_prova.py (lavoro/stazioni_prova_r8/<area>.json): stesse stazioni, stesse ore della
cassaforte (ogni 3 ore), stessa esposizione per ottante. Per ogni ora calcola con Riferimenti.kt il Fast dell'area
con il livello del radar, e sul mare davanti alla stazione (celle della libreria entro 8 km, fuori dal primo chilometro
e mezzo) il rapporto medio C / modello. La previsione è ecmwf x rapporto. Passi 3 e 6 con il profilo dell'aria (stabw),
passi 9 e 12 con l'aria mescolata, come nell'app quando il GRIB non ha il profilo.
Scrive lavoro/rev/staz_<area><SUFFISSO>.json: per riga st, t, ora, staz_mare, ecmwf, c.
"""
import datetime as dt
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cfd"))
from libreria_radar import Libreria          # noqa: E402
from riferimenti import leggi_rif, CP, MAPPE  # noqa: E402

L = "/home/claude/lavoro"; D = "/home/claude/dati"
SUF = os.environ.get("SUFFISSO", ""); PROPS = os.environ.get("PROPS", "").split(); KM = os.environ.get("KM", "")
MAX = int(os.environ.get("MAX", "400")); R_KM = 8.0
KM_AREA = {"messina": "5"}


def atlante(a):
    for p in (f"{D}/atlanti/{a}.lna", f"{D}/atlanti_nuovi/atlante_{a}.lna", f"{D}/atlanti_cfd/atlante_{a}.lna"):
        if os.path.exists(p):
            return p


def area(a):
    J = json.load(open(f"{L}/stazioni_prova_r8/{a}.json")); righe = J["righe"]
    if not righe:
        print(a, "senza confronti"); return
    lib_p = f"{L}/libc6/lra1_ecmwf_{a}.lra"; LE = Libreria(lib_p)
    la_c = LE.lat0 + (LE.celle // LE.nx) * LE.dlat; lo_c = LE.lon0 + (LE.celle % LE.nx) * LE.dlon
    staz = {s["fonte"] + "_" + s["codice"]: s for s in J["stazioni"]}
    celle = {}
    for k, s in staz.items():
        d = 111.2 * np.hypot(la_c - s["lat"], (lo_c - s["lon"]) * np.cos(np.radians(s["lat"])))
        celle[k] = np.nonzero((d <= R_KM) & (LE.dist >= 1.5))[0]
    tempi = np.array(sorted({r["t"] for r in righe}))
    if len(tempi) > MAX:
        tempi = tempi[np.random.default_rng(7).permutation(len(tempi))[:MAX]]; tempi.sort()
    scene = []
    for tv in tempi:
        d = dt.datetime.fromtimestamp(int(tv), dt.timezone.utc)
        run = d.replace(hour=0 if 3 <= d.hour <= 12 else 12, minute=0, second=0)
        if d.hour == 0:
            run = run - dt.timedelta(days=1)
        passo = int((d - run).total_seconds() // 3600)
        pw = f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2" if passo <= 6 else f"{D}/e912/e912_{run:%Y%m%d%H}.grb2"
        pp = f"{D}/stabw/stabw_{run:%Y%m%d%H}.grb2" if passo <= 6 else "nessuno"
        if os.path.exists(pw):
            scene.append((int(tv), [pw, pp, str(int(tv) * 1000), "nan", "nan"]))
    print(f"{a}{SUF}: {len(scene)} ore", flush=True)
    rapp = {}
    for k0 in range(0, len(scene), 15):
        lotto = scene[k0:k0 + 15]; usc = f"/tmp/rif/staz_{a}{SUF}_{k0}.bin"
        cmd = ["java", "-Xmx3g"] + PROPS + ["-Dlens.soloFast=true", "-cp", CP, "RiferimentiKt", atlante(a), lib_p, MAPPE, KM or KM_AREA.get(a, "30"), usc] + sum([s[1] for s in lotto], [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-1500:]); raise SystemExit(1)
        box, nx, ny, out = leggi_rif(usc); os.remove(usc)
        S, N, O, E = box
        alat = N - (np.arange(ny) + 0.5) * (N - S) / ny; alon = O + (np.arange(nx) + 0.5) * (E - O) / nx
        LA, LO = np.meshgrid(alat, alon, indexing="ij")
        Jc = np.rint((LA - LE.lat0) / LE.dlat).astype(int); Ic = np.rint((LO - LE.lon0) / LE.dlon).astype(int)
        dentro = (Jc >= 0) & (Jc < LE.ny) & (Ic >= 0) & (Ic < LE.nx)
        key = np.where(dentro, Jc * LE.nx + Ic, -1)
        for (tv, _), (t2, fr, n, F) in zip(lotto, out):
            if not np.isfinite(F[1]).any():
                continue
            ok = (F[9] > 0.5) & dentro & np.isfinite(F[1]) & (F[1] > 0.3)
            c = np.hypot(F[2], F[3]) * F[6]; a_ = np.hypot(F[2], F[3])
            for k, cel in celle.items():
                m = ok & np.isin(key, LE.celle[cel])
                if m.sum() >= 20:
                    rapp[(k, tv)] = (float(c[m].mean() / F[1][m].mean()), float(a_[m].mean() / F[1][m].mean()))
        print(f"  {min(k0 + 15, len(scene))}/{len(scene)}", flush=True)
    out = []
    for r in righe:
        q = rapp.get((r["st"], r["t"]))
        if q is not None:
            out.append(dict(st=r["st"], t=r["t"], ora=r["ora"], S0=r["S0"], staz_mare=r["staz_mare"], ecmwf=r["ecmwf"], ecmwf_c=r["ecmwf_c"], c=r["ecmwf"] * q[0], a=r["ecmwf"] * q[1]))
    json.dump(out, open(f"{L}/rev/staz_{a}{SUF}.json", "w"))
    e = lambda k: np.sqrt(np.mean([(x[k] - x["staz_mare"]) ** 2 for x in out])) / 0.5144   # noqa: E731
    b = lambda k: np.mean([(x[k] - x["staz_mare"]) for x in out]) / 0.5144                 # noqa: E731
    if out:
        print(f"{a}{SUF}: {len(out)} confronti; rms/medio nodi  modello {e('ecmwf'):.2f}/{b('ecmwf'):+.2f}  Lens2 {e('ecmwf_c'):.2f}/{b('ecmwf_c'):+.2f}  C {e('c'):.2f}/{b('c'):+.2f}  senza livello {e('a'):.2f}/{b('a'):+.2f}", flush=True)


if __name__ == "__main__":
    for a in sys.argv[1:]:
        area(a)
