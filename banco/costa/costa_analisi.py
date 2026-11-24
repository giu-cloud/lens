"""Atlante del vento sotto costa dal radar dei Sentinel-1 e confronto con ECMWF e con Lens.

  python3 costa_analisi.py <area> <cartella dei .npz> <atlante .lna o -> <cartella d'uscita>

Per ogni scena: il vento del radar (SAR, 1 km), quello del modello alla stessa ora e nello stesso punto (ECMWF, nel
prodotto), e quello di Lens fatto passare dall'atlante dell'area (LensSar.kt). Confronto sulle celle sotto costa, dopo
aver tolto a ciascun campo la sua media al largo: si misura la forma (ombre, accelerazioni), non l'errore d'insieme.
"""
import glob
import json
import os
import subprocess
import sys

import numpy as np
from scipy import ndimage

S = os.path.dirname(os.path.abspath(__file__))
KN = 1.943844


def nome(d, *cand, contiene=None, escludi=None):
    for c in cand:
        if c in d:
            return c
    if contiene:
        for k in d:
            if all(x in k for x in contiene) and not (escludi and any(e in k for e in escludi)):
                return k
    return None


def uv(d, pref=""):
    """u, v (verso dove va il vento, m/s) dalle variabili del file."""
    if pref == "":
        u = nome(d, "eastward_wind", contiene=["eastward"], escludi=["model"])
        v = nome(d, "northward_wind", contiene=["northward"], escludi=["model"])
        s = nome(d, "wind_speed", contiene=["speed"], escludi=["model"])
        di = nome(d, "wind_to_dir", "wind_dir", contiene=["dir"], escludi=["model"])
    else:
        u = nome(d, "model_eastward_wind", contiene=["model", "eastward"])
        v = nome(d, "model_northward_wind", contiene=["model", "northward"])
        s = nome(d, "model_speed", contiene=["model", "speed"])
        di = nome(d, "model_to_dir", "model_dir", contiene=["model", "dir"])
    if u and v:
        return d[u].astype("float64"), d[v].astype("float64"), (u, v)
    if s and di:
        a = np.radians(d[di].astype("float64"))           # direzione verso cui va (to_dir)
        sp = d[s].astype("float64")
        return sp * np.sin(a), sp * np.cos(a), (s, di)
    raise KeyError(f"vento {pref or 'radar'} non trovato in {list(d)}")


def carica(cartella, area):
    files = sorted(glob.glob(os.path.join(cartella, f"costa_{area}_*.npz")))
    scene = []
    lat = lon = None
    info = []
    for f in files:
        z = np.load(f, allow_pickle=False)
        d = {k: z[k] for k in z.files}
        meta = json.loads(str(d.get("meta", "{}")))
        if lat is None:
            lat, lon = d["lat"], d["lon"]
        elif d["lat"].shape != lat.shape or d["lon"].shape != lon.shape:
            print(f"{f}: griglia diversa, saltato"); continue
        su, sv, ks = uv(d)
        mu, mv, km = uv(d, "model")
        mt = d.get(nome(d, "measurement_time", contiene=["time"], escludi=["model"]) or "", None)
        q = {k: d[k] for k in d if any(x in k for x in ("quality", "flag", "info"))}
        info.append(f"{os.path.basename(f)}: {len(d['time'])} scene; radar {ks}, modello {km}; altre {sorted(set(d) - {'lat', 'lon', 'time', 'meta'} - set(ks) - set(km))}")
        for i in range(len(d["time"])):
            t = d["time"][i]
            if mt is not None and np.isfinite(mt[i]).any():
                t = int(np.nanmean(mt[i]) * 3600)
            scene.append(dict(t=int(t), su=su[i], sv=sv[i], mu=mu[i], mv=mv[i], q={k: v[i] for k, v in q.items()},
                              ds=meta.get("dataset", os.path.basename(f))))
    # stessa scena in due dataset (NRT e MY): si tiene la prima
    scene.sort(key=lambda s: s["t"])
    tenute = []
    for s in scene:
        if tenute and abs(s["t"] - tenute[-1]["t"]) < 1800:
            continue
        tenute.append(s)
    return lat, lon, tenute, info


def lens(scene, lat, lon, atlante, area, out):
    if atlante == "-" or not os.path.exists(atlante):
        return None
    ny, nx = len(lat), len(lon)
    ib = os.path.join(out, f"sar_{area}.bin"); ob = os.path.join(out, f"lens_{area}.bin")
    with open(ib, "wb") as f:
        f.write(b"SAR1")
        f.write(np.array([ny, nx], "<i4").tobytes())
        f.write(np.array([lat[0], lat[1] - lat[0], lon[0], lon[1] - lon[0]], "<f8").tobytes())
        f.write(np.array([len(scene)], "<i4").tobytes())
        for s in scene:
            f.write(np.array([s["t"] * 1000], "<i8").tobytes())
            f.write(s["mu"].astype("<f4").tobytes()); f.write(s["mv"].astype("<f4").tobytes())
    cp = f"{S}/lenssar.jar:{S}/../ut:/tmp/kotlinc/lib/kotlin-stdlib.jar"
    props = [f"-D{k}={v}" for k, v in json.loads(os.environ.get("LENS_PROPS", "{}")).items()]
    subprocess.run(["java", "-Xmx3g", *props, "-cp", cp, "LensSarKt", atlante, ib, ob], check=True)
    raw = np.fromfile(ob, dtype="<f4", offset=16)
    raw = raw.reshape(len(scene), 2, ny, nx)
    return raw[:, 0].astype("float64"), raw[:, 1].astype("float64")


def main():
    area, cartella, atlante, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    lat, lon, scene, info = carica(cartella, area)
    print("\n".join(info))
    print(f"{area}: {len(scene)} scene")
    if not scene:
        return
    ny, nx = len(lat), len(lon)
    SU = np.array([s["su"] for s in scene]); SV = np.array([s["sv"] for s in scene])
    MU = np.array([s["mu"] for s in scene]); MV = np.array([s["mv"] for s in scene])
    S_ = np.hypot(SU, SV); M_ = np.hypot(MU, MV)
    # terra: celle mai misurate dal radar; distanza dalla terra in km (celle di 0,01 gradi)
    mai = ~np.isfinite(S_).any(axis=0)
    kmy = 111.19 * abs(lat[1] - lat[0]); kmx = 111.19 * np.cos(np.radians(lat.mean())) * abs(lon[1] - lon[0])
    dist = ndimage.distance_transform_edt(~mai, sampling=(kmy, kmx))
    costa = (dist >= 1.5) & (dist <= 10.0)
    largo = dist > 10.0
    if largo.sum() < 50:
        largo = dist > np.percentile(dist[~mai], 75)
    L = lens(scene, lat, lon, atlante, area, out)
    LS = np.hypot(*L) if L else None
    # per ogni scena: direzione e forza del modello al largo, e le forme (campo / media al largo)
    righe = []
    pat = {"radar": [], "ecmwf": [], "lens": []}
    for k, s in enumerate(scene):
        ok = np.isfinite(S_[k]) & np.isfinite(M_[k])
        ol = ok & largo
        if ol.sum() < 30 or (ok & costa).sum() < 30:
            continue
        mu = np.nanmean(MU[k][ol]); mv = np.nanmean(MV[k][ol]); ms = np.nanmean(M_[k][ol])
        if ms < 3.0:
            continue
        dir_da = (np.degrees(np.arctan2(-mu, -mv)) + 360) % 360      # da dove viene
        rs = S_[k] / np.nanmean(S_[k][ol]); rm = M_[k] / ms
        rl = None
        if LS is not None and np.isfinite(LS[k][ol]).sum() > 30:
            rl = LS[k] / np.nanmean(LS[k][ol])
        c = ok & costa
        r = dict(t=s["t"], dir=float(dir_da), kn=float(ms * KN), radar_largo_kn=float(np.nanmean(S_[k][ol]) * KN),
                 err_ecmwf=float(np.sqrt(np.nanmean((rm[c] - rs[c]) ** 2))),
                 bias_ecmwf=float(np.nanmean(rm[c] - rs[c])))
        if rl is not None:
            cl = c & np.isfinite(rl)
            r["err_lens"] = float(np.sqrt(np.nanmean((rl[cl] - rs[cl]) ** 2)))
            r["bias_lens"] = float(np.nanmean(rl[cl] - rs[cl]))
        righe.append(r)
        pat["radar"].append(rs); pat["ecmwf"].append(rm); pat["lens"].append(rl if rl is not None else np.full_like(rs, np.nan))
    print(f"{area}: {len(righe)} scene usate (vento al largo almeno 6 nodi, abbastanza celle)")
    with open(os.path.join(out, f"scene_{area}.json"), "w") as f:
        json.dump(righe, f)
    if not righe:
        return
    e = np.array([r["err_ecmwf"] for r in righe]); el = np.array([r.get("err_lens", np.nan) for r in righe])
    print(f"  forma sotto costa (1,5-10 km dalla terra), errore quadratico medio sul rapporto col largo:")
    print(f"    ECMWF {np.mean(e):.3f}   Lens {np.nanmean(el):.3f}   Lens meglio in {np.nansum(el < e)} scene su {np.isfinite(el).sum()}")
    # atlante per settori di 45 gradi
    dirs = np.array([r["dir"] for r in righe])
    sett = ((dirs + 22.5) // 45).astype(int) % 8
    nomi = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    atl = {}
    for q in range(8):
        idx = np.where(sett == q)[0]
        if len(idx) < 5:
            continue
        atl[nomi[q]] = {k: np.nanmean(np.array([pat[k][i] for i in idx]), axis=0) for k in pat}
        atl[nomi[q]]["n"] = len(idx)
        ee = e[idx]; ll = el[idx]
        print(f"    da {nomi[q]:>2}: {len(idx):3d} scene, forza al largo {np.mean([righe[i]['kn'] for i in idx]):4.1f} nodi; errore ECMWF {np.mean(ee):.3f}, Lens {np.nanmean(ll):.3f}")
    np.savez_compressed(os.path.join(out, f"atlante_{area}.npz"), lat=lat, lon=lon, dist=dist,
                        **{f"{d}_{k}": v[k] for d, v in atl.items() for k in ("radar", "ecmwf", "lens")},
                        **{f"{d}_n": np.array(v["n"]) for d, v in atl.items()})
    disegna(area, lat, lon, mai, atl, out)


def disegna(area, lat, lon, mai, atl, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    dd = list(atl)
    if not dd:
        return
    fig, ax = plt.subplots(3, len(dd), figsize=(3.2 * len(dd), 8.4), squeeze=False)
    ext = [lon[0], lon[-1], lat[0], lat[-1]]
    asp = 1 / np.cos(np.radians(lat.mean()))
    for j, d in enumerate(dd):
        for i, k in enumerate(("radar", "ecmwf", "lens")):
            a = ax[i][j]
            f = np.ma.masked_invalid(atl[d][k])
            im = a.imshow(f, origin="lower", extent=ext, cmap="RdBu_r", norm=TwoSlopeNorm(1.0, 0.3, 1.5), aspect=asp)
            a.contour(lon, lat, mai.astype(float), [0.5], colors="k", linewidths=0.6)
            a.set_xticks([]); a.set_yticks([])
            if i == 0:
                a.set_title(f"vento da {d} ({atl[d]['n']} scene)", fontsize=9)
            if j == 0:
                a.set_ylabel({"radar": "radar (misurato)", "ecmwf": "ECMWF", "lens": "Lens su ECMWF"}[k], fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.6, label="vento rispetto al largo (1 = come al largo)")
    fig.suptitle(f"{area}: vento medio rispetto al largo, per direzione del vento", fontsize=11)
    fig.savefig(os.path.join(out, f"atlante_{area}.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
