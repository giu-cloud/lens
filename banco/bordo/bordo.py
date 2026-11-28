"""Il bordo dei riquadri di Lens (28/11): i dati. Per i passaggi della cassaforte di un'area grande, su ogni cella del
radar: il vento vero, il modello della tabella (mod) e il rapporto C / modello di Lens con il livello misurato su tutta
l'area e sui quattro riquadri finti (un lato rientrato di RIENTRO km). Usa Bordo.kt.

  [RIENTRO=15] python3 bordo.py <area> [<area> ...]     -> lavoro/bordo/bordo_<area>.npz
"""
import datetime as dt, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cfd"))
from riferimenti import AREE, L, D, MAPPE   # noqa: E402

CP = "/tmp/lens-motore-059/motore.jar:/tmp/bordo/bordo.jar:/opt/kotlin2021/kotlinc/lib/kotlin-stdlib.jar"
RIENTRO = os.environ.get("RIENTRO", "15")
SETTIMANA = 7 * 86400


def leggi(p):
    b = open(p, "rb").read(); assert b[:4] == b"BRD1"
    ns, nx, ny = np.frombuffer(b, "<i4", 3, 4); box = np.frombuffer(b, "<f8", 4, 16)
    o = 48; out = []
    for _ in range(ns):
        t = int(np.frombuffer(b, "<i8", 1, o)[0]); F = np.frombuffer(b, "<f4", 8 * nx * ny, o + 8).reshape(8, ny, nx); o += 8 + 32 * nx * ny
        out.append((t, F))
    return box, nx, ny, out


def area(nome):
    c = AREE[nome]
    T = np.load(f"{L}/{c['tab']}.npz", mmap_mode="r"); sel = np.asarray(T["area"]) == c["idx"]
    g = lambda k: np.asarray(T[k][sel])   # noqa: E731
    vero = np.asarray(np.load(f"{L}/radar_vero_{c['tab']}.npy", mmap_mode="r")[sel])
    G = np.load(c["griglia"]); rlat, rlon = np.asarray(G["lat"]), np.asarray(G["lon"])
    sc, t, j, i, mod, dist, v925, d925 = g("scena"), g("t"), g("j").astype(int), g("i").astype(int), g("mod"), g("dist"), g("v925"), g("d925")
    righe = np.where(((t // SETTIMANA) % 5 == 2) & np.isfinite(vero))[0]
    scene = []
    for s in np.unique(sc[righe]):
        m = righe[sc[righe] == s]; tt = int(t[m[0]]); d = dt.datetime.fromtimestamp(tt, dt.timezone.utc)
        run = d.replace(hour=0 if d.hour < 12 else 12, minute=0, second=0)
        pw = f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2"; pp = f"{D}/stabw/stabw_{run:%Y%m%d%H}.grb2"
        if not (os.path.exists(pw) and os.path.exists(pp)):
            continue
        sp = np.nanmedian(v925[m]); dd = np.radians(np.nanmedian(d925[m]))
        u, v = (-sp * np.sin(dd), -sp * np.cos(dd)) if np.isfinite(sp) and np.isfinite(dd) else (np.nan, np.nan)
        scene.append((s, tt, [pw, pp, str(tt * 1000), f"{u:.3f}", f"{v:.3f}"]))
    print(f"{nome}: {len(scene)} passaggi", flush=True)
    R = {k: [] for k in ("scena", "t", "lat", "lon", "dist", "vero", "mod", "r_full", "r_n", "r_s", "r_o", "r_e")}
    for k0 in range(0, len(scene), 15):
        lotto = scene[k0:k0 + 15]; usc = f"/tmp/bordo/{nome}_{k0}.bin"
        r = subprocess.run(["java", "-Xmx4g", "-cp", CP, "BordoKt", c["atl"], c["lib"], MAPPE, str(c["km"] if c["km"] != 30 else 10), RIENTRO, usc] + sum([s[2] for s in lotto], []), capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-2000:]); raise SystemExit(1)
        box, nx, ny, out = leggi(usc); os.remove(usc)
        S, N, O, E = box
        alat = N - (np.arange(ny) + 0.5) * (N - S) / ny; alon = O + (np.arange(nx) + 0.5) * (E - O) / nx
        LA, LO = np.meshgrid(alat, alon, indexing="ij")
        J = np.rint((LA - rlat[0]) / (rlat[1] - rlat[0])).astype(int); I = np.rint((LO - rlon[0]) / (rlon[1] - rlon[0])).astype(int)
        for (s, tt, _), (t2, F) in zip(lotto, out):
            M10, sF, mare = F[0], F[1], F[2] > 0.5
            ok = mare & (J >= 0) & (J < len(rlat)) & (I >= 0) & (I < len(rlon)) & np.isfinite(M10) & (M10 > 0.3) & np.isfinite(sF)
            key = (J * len(rlon) + I)[ok]; nb = np.bincount(key, minlength=len(rlat) * len(rlon))
            med = lambda x: np.bincount(key, x[ok], minlength=len(rlat) * len(rlon)) / np.maximum(nb, 1)   # noqa: E731
            m10 = np.maximum(med(M10), 1e-6)
            rr = [med(sF * F[3 + q]) / m10 for q in range(5)]
            m = righe[sc[righe] == s]; q = j[m] * len(rlon) + i[m]; buoni = nb[q] >= 3; m = m[buoni]; q = q[buoni]
            for kk, vv in (("scena", sc[m]), ("t", t[m]), ("lat", rlat[j[m]]), ("lon", rlon[i[m]]), ("dist", dist[m]), ("vero", vero[m]), ("mod", mod[m]),
                           ("r_full", rr[0][q]), ("r_n", rr[1][q]), ("r_s", rr[2][q]), ("r_o", rr[3][q]), ("r_e", rr[4][q])):
                R[kk].append(np.asarray(vv))
        print(f"  {min(k0 + 15, len(scene))}/{len(scene)}", flush=True)
    os.makedirs(f"{L}/bordo", exist_ok=True)
    np.savez_compressed(f"{L}/bordo/bordo_{nome}.npz", box=np.array(c.get("box_area", [])), **{k: np.concatenate(v) for k, v in R.items()})


if __name__ == "__main__":
    for a in sys.argv[1:]:
        area(a)
