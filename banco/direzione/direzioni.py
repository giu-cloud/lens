"""La direzione di Lens contro le stazioni sulla costa, a tutte le ore (27/11 sera). Usa Direzioni.kt.

  [N=700] [ATL=<atlante>] [SUF=_fine] python3 direzioni.py <area> [<area> ...]

Stazioni dell'elenco (lavoro/stazioni_elenco.json) sotto i 200 m. Ore: ogni 3 ore (passi 3, 6, 9, 12 di ECMWF), N a caso
tra tutte le settimane; le settimane della cassaforte servono per la verifica, le altre per tarare.
Scrive lavoro/dir/dir_<area><SUF>.npz: per riga stazione, tempo, direzione e forza della stazione, e u, v di M, L0, LA, LM, LS
nella cella della stazione (c) e nella cella di mare aperto più vicina (m).
"""
import datetime as dt, json, os, subprocess, sys
import numpy as np

L = "/home/claude/lavoro"; D = "/home/claude/dati"
CP = "/tmp/lens-motore-059/motore.jar:/tmp/dir/dir.jar:/opt/kotlin2021/kotlinc/lib/kotlin-stdlib.jar"
MAPPE = "/home/claude/shift/app/src/main/assets/maps"
N = int(os.environ.get("N", "700")); SUF = os.environ.get("SUF", "")


def atlante(a):
    if os.environ.get("ATL"):
        return os.environ["ATL"]
    for p in (f"{D}/atlanti/{a}.lna", f"{D}/atlanti_nuovi/atlante_{a}.lna", f"{D}/atlanti_cfd/atlante_{a}.lna"):
        if os.path.exists(p):
            return p


def file_staz(s):
    for d in ("stazioni_r4", "stazioni_r3"):
        p = f"{D}/{d}/{s['fonte']}_{s['codice']}.npz"
        if os.path.exists(p):
            return p


def area(a):
    E = [s for s in json.load(open(f"{L}/stazioni_elenco.json")) if s["area"] == a and not (s.get("quota") or 0) > 200 and file_staz(s)]
    if not E or not atlante(a):
        print(a, "senza stazioni o senza atlante"); return
    S = {}
    for s in E:
        Z = np.load(file_staz(s)); o = np.argsort(Z["t"]); S[s["fonte"] + "_" + s["codice"]] = (Z["t"][o], Z["vento"][o], Z["direzione"][o])
    t0 = min(x[0][0] for x in S.values()); t1 = max(x[0][-1] for x in S.values())
    tt = np.arange(int(np.ceil(t0 / 10800)) * 10800, t1, 10800)
    tt = tt[np.random.default_rng(11).permutation(len(tt))]
    scene = []
    for tv in tt:
        d = dt.datetime.fromtimestamp(int(tv), dt.timezone.utc)
        run = d.replace(hour=0 if 3 <= d.hour <= 12 else 12, minute=0, second=0)
        if d.hour == 0:
            run = run - dt.timedelta(days=1)
        passo = int((d - run).total_seconds() // 3600)
        pw = f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2" if passo <= 6 else f"{D}/e912/e912_{run:%Y%m%d%H}.grb2"
        pp = f"{D}/stabw/stabw_{run:%Y%m%d%H}.grb2" if passo <= 6 else "nessuno"
        if os.path.exists(pw):
            scene.append((pw, pp, int(tv)))
        if len(scene) >= N:
            break
    scene.sort()
    punti = ",".join(f"{k}:{s['lat']}:{s['lon']}" for k, s in zip(S, E))
    righe = []
    os.makedirs(f"{L}/dir", exist_ok=True)
    for k0 in range(0, len(scene), 80):
        usc = f"/tmp/dir/{a}{SUF}_{k0}.txt"
        cmd = ["java", "-Xmx3g", "-cp", CP, "DirezioniKt", atlante(a), MAPPE, usc, punti] + sum([[x[0], x[1], str(x[2] * 1000)] for x in scene[k0:k0 + 80]], [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-1500:]); raise SystemExit(1)
        for ln in open(usc):
            if ln.startswith("#"):
                intest = ln.strip(); continue
            f = ln.split(); tv = int(f[0]) // 1000; ts, vs, ds = S[f[1]]
            j = np.searchsorted(ts, tv); best = None
            for q in (j - 1, j):
                if 0 <= q < len(ts) and abs(ts[q] - tv) <= 1200 and np.isfinite(vs[q]) and np.isfinite(ds[q]) and (best is None or abs(ts[q] - tv) < abs(ts[best] - tv)):
                    best = q
            if best is not None:
                righe.append([f[1], tv, vs[best], ds[best]] + [float(x.replace(",", ".")) for x in f[2:]])
        os.remove(usc)
        print(f"  {a}{SUF} {min(k0 + 80, len(scene))}/{len(scene)}", flush=True)
    st = np.array([r[0] for r in righe]); X = np.array([r[1:] for r in righe], float)
    np.savez_compressed(f"{L}/dir/dir_{a}{SUF}.npz", st=st, X=X, intest=intest)
    print(a, SUF, len(righe), "righe", intest, flush=True)


if __name__ == "__main__":
    for a in sys.argv[1:]:
        area(a)
