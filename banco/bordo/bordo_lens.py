"""La prova del bordo con il codice del motore (BordoLens.kt): modello, LensWind di prima, LensWind con l'alone, contro il
radar, sui riquadri finti (un lato rientrato di 15 km), passaggi della cassaforte.
  python3 bordo_lens.py <area> [...]      -> lavoro/bordo/brl_<area>.npz, poi il riassunto con  python3 bordo_lens.py --riassunto
"""
import os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cfd"))
from riferimenti import AREE, L, MAPPE      # noqa: E402
from bordo import SETTIMANA                  # noqa: E402
from analisi import BOX, lato               # noqa: E402
import bordo as B                            # noqa: E402

CP = "/tmp/lens-motore-0592/motore.jar:/tmp/bordo/bordolens.jar:/opt/kotlin2021/kotlinc/lib/kotlin-stdlib.jar"
KN = 1 / 0.5144


def area(nome):
    c = AREE[nome]
    T = np.load(f"{L}/{c['tab']}.npz", mmap_mode="r"); sel = np.asarray(T["area"]) == c["idx"]
    g = lambda k: np.asarray(T[k][sel])   # noqa: E731
    vero = np.asarray(np.load(f"{L}/radar_vero_{c['tab']}.npy", mmap_mode="r")[sel])
    G = np.load(c["griglia"]); rlat, rlon = np.asarray(G["lat"]), np.asarray(G["lon"])
    sc, t, j, i = g("scena"), g("t"), g("j").astype(int), g("i").astype(int)
    # le stesse scene di bordo.py
    Z = np.load(f"{L}/bordo/bordo_{nome}.npz"); tt = np.unique(Z["t"])
    righe = np.where(np.isin(t, tt) & np.isfinite(vero))[0]
    import datetime as dt
    D = "/home/claude/dati"; scene = []
    for s in np.unique(sc[righe]):
        m = righe[sc[righe] == s]; t0 = int(t[m[0]]); d = dt.datetime.fromtimestamp(t0, dt.timezone.utc)
        run = d.replace(hour=0 if d.hour < 12 else 12, minute=0, second=0)
        v925, d925 = g("v925")[m], g("d925")[m]
        sp = np.nanmedian(v925); dd = np.radians(np.nanmedian(d925))
        u, v = (-sp * np.sin(dd), -sp * np.cos(dd)) if np.isfinite(sp) and np.isfinite(dd) else (np.nan, np.nan)
        scene.append((s, t0, [f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2", f"{D}/stabw/stabw_{run:%Y%m%d%H}.grb2", str(t0 * 1000), f"{u:.3f}", f"{v:.3f}"]))
    out = {k: [] for k in ("t", "lato", "lat", "lon", "dist", "vero", "M", "prima", "alone")}
    for k0 in range(0, len(scene), 15):
        lotto = scene[k0:k0 + 15]; usc = f"/tmp/bordo/brl_{nome}_{k0}.bin"
        r = subprocess.run(["java", "-Xmx4g", "-cp", CP, "BordoLensKt", c["atl"], c["lib"], MAPPE, str(c["km"] if c["km"] != 30 else 10), "15", usc] + sum([s[2] for s in lotto], []), capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-2000:]); raise SystemExit(1)
        b = open(usc, "rb").read(); os.remove(usc); assert b[:4] == b"BRL1"
        ns, ny, nx = np.frombuffer(b, "<i4", 3, 4); lat0, lon0, dla, dlo = np.frombuffer(b, "<f8", 4, 16)
        assert abs(lat0 - rlat[0]) < 1e-3 and abs(lon0 - rlon[0]) < 1e-3 and ny == len(rlat) and nx == len(rlon), (lat0, rlat[0], ny, len(rlat))
        o = 48
        for (s, t0, _) in lotto:
            o += 8; F = np.frombuffer(b, "<f4", 12 * ny * nx, o).reshape(4, 3, ny, nx); o += 4 * 12 * ny * nx
            m = righe[sc[righe] == s]
            for q, nm in enumerate("nsoe"):
                for kk, vv in (("t", t[m]), ("lato", np.full(len(m), nm)), ("lat", rlat[j[m]]), ("lon", rlon[i[m]]), ("dist", g("dist")[m]), ("vero", vero[m]),
                               ("M", F[q, 0, j[m], i[m]]), ("prima", F[q, 1, j[m], i[m]]), ("alone", F[q, 2, j[m], i[m]])):
                    out[kk].append(np.asarray(vv))
        print(f"  {nome} {min(k0 + 15, len(scene))}/{len(scene)}", flush=True)
    np.savez_compressed(f"{L}/bordo/brl_{nome}.npz", **{k: np.concatenate(v) for k, v in out.items()})


def riassunto():
    rng = np.random.default_rng(4); R = []
    for a in BOX:
        p = f"{L}/bordo/brl_{a}.npz"
        if not os.path.exists(p):
            continue
        Z = np.load(p); X = {k: Z[k] for k in Z.files}
        for s in "nsoe":
            k = X["lato"] == s
            Y = {kk: v[k] for kk, v in X.items()}
            d, x, altri = lato(Y, BOX[a], s)
            ok = np.isfinite(Y["vero"]) & (Y["vero"] > 0.5) & (Y["M"] > 1.0) & np.isfinite(Y["prima"]) & np.isfinite(Y["alone"]) & (d > -15) & (d < 8) & (altri > 3)
            R.append(dict(area=a, d=d[ok], t=Y["t"][ok], **{kk: Y[kk][ok] for kk in ("vero", "M", "prima", "alone")}))
    print("Errore medio / quadratico (nodi) per fascia di distanza dal lato finto, con il codice del motore:")
    for lo_, hi in [(-15, -10), (-10, -5), (-5, -2), (-2, 0), (0, 2.8), (2.8, 8)]:
        e = {n: np.concatenate([(r[n] - r["vero"])[(r["d"] >= lo_) & (r["d"] < hi)] for r in R]) * KN for n in ("M", "prima", "alone")}
        print(f"  {lo_:5.1f}..{hi:4.1f} {len(e['M']):6d}  " + "  ".join(f"{n}: {e[n].mean():+5.2f}/{np.sqrt((e[n] ** 2).mean()):4.2f}" for n in e))
    print("Per area, da 10 km fuori a 2,8 km dentro: quadratico prima, con l'alone, differenza (intervallo 90%)")
    for a in list(dict.fromkeys(r["area"] for r in R)) + ["TUTTE"]:
        rr = [r for r in R if a == "TUTTE" or r["area"] == a]
        sel = [(r, (r["d"] >= -10) & (r["d"] < 2.8)) for r in rr]
        e1 = np.concatenate([((r["prima"] - r["vero"]) * KN)[m] for r, m in sel]); e2 = np.concatenate([((r["alone"] - r["vero"]) * KN)[m] for r, m in sel])
        tt = np.concatenate([r["t"][m] for r, m in sel]); U = np.unique(tt); pos = {u: np.where(tt == u)[0] for u in U}; bb = []
        for _ in range(300):
            p = np.concatenate([pos[u] for u in rng.choice(U, len(U))])
            bb.append(np.sqrt((e2[p] ** 2).mean()) - np.sqrt((e1[p] ** 2).mean()))
        q1, q2 = np.sqrt((e1 ** 2).mean()), np.sqrt((e2 ** 2).mean())
        print(f"  {a:13s} {len(e1):6d}  {q1:.2f}  {q2:.2f}  {q2 - q1:+.3f} ({np.percentile(bb, 5):+.3f} … {np.percentile(bb, 95):+.3f})")


if __name__ == "__main__":
    if sys.argv[1] == "--riassunto":
        riassunto()
    else:
        for a in sys.argv[1:]:
            area(a)
