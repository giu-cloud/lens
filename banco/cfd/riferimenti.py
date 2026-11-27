"""I riferimenti della prova CFD (27/11, lens/cfd_piano_26_11.md, fase 0): M, Fast, Detailed e B sui passaggi del radar
delle settimane della cassaforte, sette aree, con la stabilità vera della scena (profilo ECMWF). Usa Riferimenti.kt.

  python3 riferimenti.py <area> [<area> ...]      (aree: palermo messina eolie monaco corsica-nord dalmazia cagliari)

Scrive /home/claude/lavoro/cfd_rif/rif_<area>.npz, una riga per cella del radar e passaggio (solo dentro il riquadro
dell'area fine della prova), con il vento vero del radar, il modello della tabella (mod) e i rapporti di ogni versione sul
modello lisciato che Lens riceve: le previsioni sono mod x rapporto, come nel banco di misura (stretto.py). Stampa gli
errori come lo studio dello Stretto: quadratico / medio / sistematico per cella, tutto ed entro 3 km.
"""
import datetime as dt
import os
import subprocess
import sys

import numpy as np

L = "/home/claude/lavoro"; D = "/home/claude/dati"
SETTIMANA = 7 * 86400
AREE = {
    "palermo": dict(tab="tabella6", idx=2, griglia=f"{L}/v_scia/campi_palermo.npz", atl=f"{D}/atlanti/palermo.lna", lib=f"{L}/libc6/lra1_ecmwf_palermo.lra", km=30, box=(38.0, 38.287, 13.108, 13.612)),
    "messina": dict(tab="tabella6", idx=3, griglia=f"{L}/v_scia/campi_messina.npz", atl=f"{D}/atlanti/messina.lna", lib=f"{L}/libc6/lra1_ecmwf_messina.lra", km=5, box=None),
    "eolie": dict(tab="tabella6", idx=1, griglia=f"{L}/v_scia/campi_eolie.npz", atl=f"{D}/atlanti/eolie.lna", lib=f"{L}/libc6/lra1_ecmwf_eolie.lra", km=30, box=(38.29, 38.7, 14.72, 15.12)),
    "monaco": dict(tab="tabella6", idx=12, griglia=f"{D}/sar/costa_monaco_nrt_l3-s1a-sar-asc.npz", atl=f"{D}/atlanti_cfd/atlante_monaco.lna", lib=f"{L}/libc6/lra1_ecmwf_monaco.lra", km=30, box=(43.55, 43.8, 7.15, 7.5)),
    "corsica-nord": dict(tab="tabella6", idx=11, griglia=f"{D}/sar/costa_corsica-nord_nrt_l3-s1a-sar-asc.npz", atl=f"{D}/atlanti_cfd/atlante_corsica-nord.lna", lib=f"{L}/libc6/lra1_ecmwf_corsica-nord.lra", km=30, box=(42.75, 43.1, 9.18, 9.62)),
    "dalmazia": dict(tab="tabella_n2", idx=17, griglia=f"{D}/sar/costa_dalmazia_nrt_l3-s1a-sar-asc.npz", atl=f"{D}/atlanti_cfd/atlante_dalmazia.lna", lib=f"{L}/libc_n2/lra1_ecmwf_dalmazia.lra", km=30, box=(43.25, 43.55, 16.1, 16.7)),
    "cagliari": dict(tab="tabella6", idx=6, griglia=f"{L}/v_scia/campi_cagliari.npz", atl=f"{D}/atlanti/cagliari.lna", lib=f"{L}/libc6/lra1_ecmwf_cagliari.lra", km=30, box=(38.95, 39.25, 9.0, 9.42)),
}
CP = "/tmp/lens-motore-build/motore.jar:/tmp/rif/rif.jar:/opt/kotlin2021/kotlinc/lib/kotlin-stdlib.jar"
MAPPE = "/home/claude/shift/app/src/main/assets/maps"
RIL = os.environ.get("RILIEVO", "")          # la manopola del rilievo (vuoto = quella del motore, 1)
SUF = f"_r{RIL}" if RIL else ""


def leggi_rif(p):
    b = open(p, "rb").read()
    assert b[:4] == b"RIF1"
    ns, nx, ny = np.frombuffer(b, "<i4", 3, 4); box = np.frombuffer(b, "<f8", 4, 16)
    o = 48; scene = []
    for _ in range(ns):
        t = int(np.frombuffer(b, "<i8", 1, o)[0]); fr, n = np.frombuffer(b, "<f8", 2, o + 8)
        F = np.frombuffer(b, "<f4", 12 * nx * ny, o + 24).reshape(12, ny, nx); o += 24 + 4 * 12 * nx * ny
        scene.append((t, fr, n, F))
    return box, nx, ny, scene


def sist(cel, e, mk):
    ids, inv = np.unique(cel[mk], return_inverse=True); n = np.bincount(inv); me = np.bincount(inv, e[mk]) / n
    return float(np.sqrt(np.sum(n * me ** 2) / n.sum()))


def area(nome):
    c = AREE[nome]
    T = np.load(f"{L}/{c['tab']}.npz", mmap_mode="r"); sel = np.asarray(T["area"]) == c["idx"]
    g = lambda k: np.asarray(T[k][sel])                                   # noqa: E731
    vero = np.load(f"{L}/radar_vero_{c['tab']}.npy", mmap_mode="r")[sel]
    G = np.load(c["griglia"]); rlat, rlon = np.asarray(G["lat"]), np.asarray(G["lon"])
    sc, t, j, i, mod, dist, v925, d925 = g("scena"), g("t"), g("j").astype(int), g("i").astype(int), g("mod"), g("dist"), g("v925"), g("d925")
    la, lo = rlat[j], rlon[i]
    cass = (t // SETTIMANA) % 5 == 2
    dentro = np.ones(len(t), bool) if c["box"] is None else (la >= c["box"][0]) & (la <= c["box"][1]) & (lo >= c["box"][2]) & (lo <= c["box"][3])
    righe = np.where(cass & dentro & np.isfinite(vero))[0]
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
    print(f"{nome}{SUF}: {len(scene)} passaggi della cassaforte con GRIB e profilo", flush=True)
    os.makedirs(f"{L}/cfd_rif", exist_ok=True)
    righe_out = {k: [] for k in ("scena", "t", "j", "i", "lat", "lon", "dist", "vero", "mod", "Fr", "rFast", "rDet", "rB", "rAf", "rAd", "dirFast", "dirDet", "dirM")}
    for k0 in range(0, len(scene), 15):
        lotto = scene[k0:k0 + 15]
        usc = f"/tmp/rif/{nome}{SUF}_{k0}.bin"
        cmd = ["java", "-Xmx5g"] + (["-Dshift.lens.relief=" + RIL, "-Dlens.soloFast=true"] if RIL else []) + ["-cp", CP, "RiferimentiKt", c["atl"], c["lib"], MAPPE, str(c["km"]), usc] + sum([s[2] for s in lotto], [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-2000:]); raise SystemExit(1)
        box, nx, ny, out = leggi_rif(usc)
        os.remove(usc)
        S, N, O, E = box
        alat = N - (np.arange(ny) + 0.5) * (N - S) / ny; alon = O + (np.arange(nx) + 0.5) * (E - O) / nx
        LA, LO = np.meshgrid(alat, alon, indexing="ij")
        dj = rlat[1] - rlat[0]; di = rlon[1] - rlon[0]
        J = np.rint((LA - rlat[0]) / dj).astype(int); I = np.rint((LO - rlon[0]) / di).astype(int)
        for (s, tt, _), (t2, fr, n, F) in zip(lotto, out):
            mare = F[9] > 0.5
            M10 = F[1]
            sF = np.hypot(F[2], F[3]) * F[6]; sD = np.hypot(F[4], F[5]) * F[7]; sB = np.hypot(F[2], F[3]) * F[8]
            sAf = np.hypot(F[2], F[3]); sAd = np.hypot(F[4], F[5])
            ok = mare & (J >= 0) & (J < len(rlat)) & (I >= 0) & (I < len(rlon)) & np.isfinite(M10) & (M10 > 0.3)
            key = (J * len(rlon) + I)[ok]; nb = np.bincount(key, minlength=len(rlat) * len(rlon))
            def media(x):
                return np.bincount(key, x[ok], minlength=len(rlat) * len(rlon)) / np.maximum(nb, 1)
            def media_dir(u, v):
                return np.degrees(np.arctan2(-media(u), -media(v))) % 360
            m10 = media(M10)
            rF, rD, rB, rAf, rAd = (media(x) / np.maximum(m10, 1e-6) for x in (sF, sD, sB, sAf, sAd))
            dF = media_dir(F[2], F[3]); dD = media_dir(F[4], F[5])
            m = righe[sc[righe] == s]
            q = j[m] * len(rlon) + i[m]; buoni = nb[q] >= 3
            m = m[buoni]; q = q[buoni]
            for kk, vv in (("scena", sc[m]), ("t", t[m]), ("j", j[m]), ("i", i[m]), ("lat", la[m]), ("lon", lo[m]), ("dist", dist[m]),
                           ("vero", vero[m]), ("mod", mod[m]), ("Fr", np.full(len(m), fr)), ("rFast", rF[q]), ("rDet", rD[q]), ("rB", rB[q]),
                           ("rAf", rAf[q]), ("rAd", rAd[q]), ("dirFast", dF[q]), ("dirDet", dD[q]), ("dirM", g("dir")[m])):
                righe_out[kk].append(np.asarray(vv))
        print(f"  {min(k0 + 15, len(scene))}/{len(scene)}", flush=True)
    R = {k: np.concatenate(v) if v else np.zeros(0) for k, v in righe_out.items()}
    np.savez_compressed(f"{L}/cfd_rif/rif_{nome}{SUF}.npz", **R)
    # gli errori, come lo studio dello Stretto
    vero_, mod_ = R["vero"], R["mod"]
    okk = np.isfinite(vero_) & (vero_ > 0.5) & (mod_ > 2.0) & np.isfinite(R["rFast"]) & np.isfinite(R["rDet"]) & np.isfinite(R["rB"])
    cel = R["j"] * 1000 + R["i"]
    prev = {"M": mod_, "Fast": mod_ * R["rFast"], "Detailed": mod_ * R["rDet"], "B": mod_ * R["rB"]}
    for fascia, f in (("tutto", okk), ("entro 3 km", okk & (R["dist"] < 3))):
        print(f"  {fascia:10s} " + "  ".join(f"{k}: {np.sqrt(np.mean(((p - vero_) / 0.5144)[f] ** 2)):.2f}/{np.mean(((p - vero_) / 0.5144)[f]):+.2f}/{sist(cel, (p - vero_) / 0.5144, f):.2f}"
                                             for k, p in prev.items()) + f"  ({f.sum()} righe, {len(np.unique(R['scena'][f]))} passaggi)")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        area(a)
