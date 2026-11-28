"""Dalle immagini ridotte (release strisce-<area>) alle direzioni del vento per finestra, con il modello accanto.

  python3 analizza.py <area>        -> lavoro/strisce/finestre_<area>.npz

Per ogni finestra buona: tempo, lat, lon, distanza dalla costa (km), qualità, direzione del radar (da dove viene il vento,
con l'ambiguità di 180 gradi sciolta con ECMWF), direzione e forza di ECMWF (GRIB a 0,25 gradi, passi 3 e 6) e di ICON-2I
(quando c'è), e l'angolo delle strisce rispetto agli assi dell'immagine (per riconoscere le cuciture).
"""
import datetime as dt, glob, json, os, sys
import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strisce import direzioni

D = "/home/claude/dati"; L = "/home/claude/lavoro/strisce"
TERRA = {"corsica-nord": (f"{D}/strisce/terra_corsica-nord.bin", 42.35, 8.15, 0.0005)}
FASCIA_KM = 1.0


def terra(area):
    p, s, o, passo = TERRA[area]
    raw = open(p, "rb").read(); ny, nx = np.frombuffer(raw, ">i4", 2)
    m = np.frombuffer(raw, np.uint8, offset=8).reshape(ny, nx).astype(bool)
    lat = s + (np.arange(ny) + 0.5) * passo
    km = np.array([passo * 111.19, passo * 111.19 * np.cos(np.radians(lat.mean()))])
    dist = ndi.distance_transform_edt(~m, sampling=km)      # km dal punto di terra più vicino
    return m, dist, s, o, passo


def su_griglia(G, s, o, passo, lat, lon, fuori):
    j = np.floor((lat - s) / passo).astype(int); i = np.floor((lon - o) / passo).astype(int)
    ok = (j >= 0) & (j < G.shape[0]) & (i >= 0) & (i < G.shape[1])
    out = np.full(lat.shape, fuori, dtype=G.dtype if G.dtype != bool else bool)
    out[ok] = G[j[ok], i[ok]]
    return out


_grib = {}


def ecmwf(t, lat, lon):
    """u, v a 10 m di ECMWF al tempo t (unix) nei punti, dal GRIB della corsa 00 o 12 (passi 3 e 6), lineare in tempo e spazio."""
    import eccodes
    d = dt.datetime.fromtimestamp(int(t), dt.timezone.utc); run = d.replace(hour=0 if d.hour < 12 else 12, minute=0, second=0, microsecond=0)
    p = f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2"
    if not os.path.exists(p):
        return None
    if p not in _grib:
        F = {}
        with open(p, "rb") as fh:
            while True:
                h = eccodes.codes_grib_new_from_file(fh)
                if h is None:
                    break
                sn = eccodes.codes_get(h, "shortName"); st = eccodes.codes_get(h, "step")
                if sn in ("10u", "10v"):
                    ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
                    la0, lo0 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees"), eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees")
                    la1, lo1 = eccodes.codes_get(h, "latitudeOfLastGridPointInDegrees"), eccodes.codes_get(h, "longitudeOfLastGridPointInDegrees")
                    v = eccodes.codes_get_values(h).reshape(nj, ni)
                    las = np.linspace(la0, la1, nj); los = np.linspace(lo0, lo1, ni)
                    if las[0] > las[-1]:
                        las = las[::-1]; v = v[::-1]
                    F[(sn, st)] = RegularGridInterpolator((las, los), v, bounds_error=False, fill_value=np.nan)
                eccodes.codes_release(h)
        _grib.clear(); _grib[p] = F
    F = _grib[p]; ore = (d - run).total_seconds() / 3600.0
    passi = sorted({k[1] for k in F})
    a = max([x for x in passi if x <= ore], default=passi[0]); b = min([x for x in passi if x >= ore], default=passi[-1])
    w = 0.0 if a == b else (ore - a) / (b - a)
    pts = np.stack([lat, lon], -1)
    return tuple((1 - w) * F[(k, a)](pts) + w * F[(k, b)](pts) for k in ("10u", "10v"))


_icon = {}


def icon(area, t, lat, lon):
    if area not in _icon:
        p = f"{D}/icon2i_3h/icon2i_{area}.npz"
        _icon[area] = np.load(p) if os.path.exists(p) else None
    Z = _icon[area]
    if Z is None:
        return None
    it = Z["t"]; k = np.searchsorted(it, t)
    if k == 0 or k >= len(it) or it[k] - it[k - 1] > 3 * 3600 + 1:
        return None
    w = (t - it[k - 1]) / (it[k] - it[k - 1]); out = []
    for nome in ("u", "v"):
        f = (1 - w) * Z[nome][k - 1].astype(float) + w * Z[nome][k].astype(float)
        out.append(RegularGridInterpolator((Z["lat"], Z["lon"]), f, bounds_error=False, fill_value=np.nan)(np.stack([lat, lon], -1)))
    return tuple(out)


def da_dove(u, v):
    return np.degrees(np.arctan2(-u, -v)) % 360


def immagine(area, f, M):
    m, dist, s, o, passo = M
    Z = np.load(f); meta = json.loads(str(Z["meta"])); t = int(os.path.basename(f).split("_")[1])
    amp = Z["amp"].astype(np.float32); amp[amp == 0] = np.nan
    ny, nx = amp.shape
    gr = np.arange(0, ny, 25); gc = np.arange(0, nx, 25)
    fla = RegularGridInterpolator((gr, gc), Z["lat"].astype(float), bounds_error=False, fill_value=None)
    flo = RegularGridInterpolator((gr, gc), Z["lon"].astype(float), bounds_error=False, fill_value=None)
    # terra e fascia costiera fuori (alla risoluzione di 40 m, a blocchi per non riempire la memoria)
    R, C = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    la = fla((R, C)); lo = flo((R, C))
    dk = su_griglia(dist, s, o, passo, la, lo, 99.0)
    amp[dk < FASCIA_KM] = np.nan
    d = direzioni(amp)
    r, c = d["r"], d["c"]
    la_w = fla((r, c)); lo_w = flo((r, c)); dk_w = su_griglia(dist, s, o, passo, la_w, lo_w, 99.0)
    # gli assi dell'immagine sulla carta: di quanto ci si sposta (est, nord, in km) per una cella lungo le colonne e lungo le righe
    kx = 111.19 * np.cos(np.radians(la_w))
    ec = np.stack([(flo((r, c + 10)) - flo((r, c - 10))) * kx, (fla((r, c + 10)) - fla((r, c - 10))) * 111.19], -1) / 20
    er = np.stack([(flo((r + 10, c)) - flo((r - 10, c))) * kx, (fla((r + 10, c)) - fla((r - 10, c))) * 111.19], -1) / 20
    a = np.radians(d["ang"])
    vec = np.cos(a)[..., None] * ec + np.sin(a)[..., None] * er          # la direzione delle strisce, in est e nord
    asse = np.degrees(np.arctan2(vec[..., 0], vec[..., 1])) % 180        # come angolo di bussola, modulo 180
    uv = ecmwf(t, la_w, lo_w)
    if uv is None:
        return None
    dm = da_dove(*uv); sm = np.hypot(*uv)
    # l'ambiguità di 180 gradi: il verso più vicino al modello
    dif = (asse - dm + 90) % 180 - 90
    dsar = (dm + dif) % 360
    ic = icon(area, t, la_w, lo_w)
    di = da_dove(*ic) if ic is not None else np.full(la_w.shape, np.nan)
    ok = (d["n"] > 0.7) & np.isfinite(dm)
    k = ok.ravel()
    col = lambda x: np.asarray(x, float).ravel()[k]
    return dict(t=np.full(k.sum(), t, float), lat=col(la_w), lon=col(lo_w), dist=col(dk_w), q=col(d["q"]), n=col(d["n"]), sar=col(dsar), dif=col(dif),
                ecm=col(dm), ecm_v=col(sm), icon=col(di), ang_img=col(d["ang"]), img=np.full(k.sum(), hash(os.path.basename(f)) % 10**6, float)), dict(amp=amp, d=d, la=la_w, lo=lo_w, asse=asse, dm=dm, dsar=dsar, dk=dk_w, meta=meta, t=t, ec=ec, er=er, di=di)


if __name__ == "__main__":
    area = sys.argv[1]; os.makedirs(L, exist_ok=True)
    M = terra(area); righe = []
    for f in sorted(glob.glob(f"{D}/strisce/{area}_*.npz")):
        r = immagine(area, f, M)
        if r is None:
            print(os.path.basename(f), "senza modello"); continue
        righe.append(r[0])
        x = r[0]; lar = x["dist"] > 20
        print(os.path.basename(f), len(x["t"]), "finestre; al largo qualità mediana %.2f, |radar - ECMWF| mediana %.1f" % (np.median(x["q"][lar]) if lar.any() else np.nan, np.median(np.abs(x["dif"][lar])) if lar.any() else np.nan), flush=True)
    out = {k: np.concatenate([r[k] for r in righe]) for k in righe[0]}
    np.savez_compressed(f"{L}/finestre_{area}.npz", **out)
