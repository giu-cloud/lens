"""Tabella delle celle per lo studio statistico di Lens: una riga per cella del radar e per passaggio.

  python3 tabella_celle.py <cartella del metro (campi_*.npz, sho_*.bin)> <rilievo di Shift .rgr> <uscita.npz>
        [cartella con Lens 0.54] [altre cartelle di varianti nome=cartella ...]

Per ogni cella e per 16 direzioni di provenienza si misura il rilievo intorno (sul rilievo di Shift a 464 m):
  sopravvento   terra (km alla prima terra, 30 = mare aperto), h15 (quota massima entro 15 km dalla costa), ang (angolo
                massimo sotto cui si vede il rilievo), hc (quota media dei primi 3 km di terra), valle (rilievo ai lati
                sopra la linea del vento, 2-6 km dentro terra: positivo = valle o passo), f5/f15 (frazione di terra lungo
                la linea entro 5 e 15 km), angw (angolo massimo nel ventaglio di ±22,5°), fw15 (terra nel ventaglio)
  sottovento    sotto (km alla prima terra davanti), hsotto (quota massima entro 10 km da lì), angsotto (angolo)
  ai lati       lsx, ldx (km alla prima terra a sinistra e a destra, guardando dove va il vento), canale = lsx + ldx
  intorno       h10 (quota massima entro 10 km), hm10 (quota media della terra entro 10 km), terra10 (frazione di terra)
Le caratteristiche si interpolano tra i due settori vicini secondo la direzione del modello letto dal mare nella cella.
Della scena: forza e direzione del modello al largo (celle a più di 8 km da terra), rapporto radar/modello al largo,
ora e mese. Il bersaglio è il radar; il modello è quello letto dal mare; Lens 0.54 e le varianti fisiche (se date)
entrano come rapporto con il modello letto dal mare.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rilievo import Rilievo

AREE_TUTTE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]
AREE = os.environ.get("AREE", " ".join(AREE_TUTTE)).split()
NSET = 16
PASSO = 0.25
RMAX = 30.0
VMIN = 3.0          # m/s: sotto, il rapporto è solo rumore
DMAX = 25.0         # km dalla terra
LARGO = 8.0         # km: oltre, la cella conta come largo della scena


def campi(cart, area):
    z = np.load(os.path.join(cart, f"campi_{area}.npz"))
    raw = np.fromfile(os.path.join(cart, f"sho_{area}.bin"), dtype="<f4", offset=16)
    ny, nx = len(z["lat"]), len(z["lon"])
    return z, raw.reshape(-1, 8, ny, nx)


def statiche(lat, lon, ril):
    s, n, w, e = lat.min() - 0.4, lat.max() + 0.4, lon.min() - 0.5, lon.max() + 0.5
    rl, ro, H = ril.griglia(0, s, n, w, e)
    H = np.maximum(H, 0).astype("f4")
    c = rl[1] - rl[0]
    kmy = 111.19; kmx = 111.19 * np.cos(np.radians(lat.mean()))
    LA, LO = np.meshgrid(lat.astype("f8"), lon.astype("f8"), indexing="ij")

    def quota(la, lo):
        j = np.clip(np.round((la - rl[0]) / c).astype(int), 0, len(rl) - 1)
        i = np.clip(np.round((lo - ro[0]) / c).astype(int), 0, len(ro) - 1)
        return H[j, i]

    passi = np.arange(PASSO, RMAX + 1e-9, PASSO)
    r3 = passi[:, None, None]

    def raggio(dy, dx):
        return np.stack([quota(LA + r * dy / kmy, LO + r * dx / kmx) for r in passi])

    def prima(hr):
        t = hr > 0
        p = np.where(t.any(axis=0), t.argmax(axis=0), len(passi))
        return np.where(p < len(passi), passi[np.minimum(p, len(passi) - 1)], RMAX)

    nomi = ["terra", "h15", "ang", "hc", "valle", "f5", "f15", "sotto", "hsotto", "angsotto", "lsx", "ldx"]
    out = {k: np.zeros((NSET,) + LA.shape, dtype="f4") for k in nomi}
    for k in range(NSET):
        a = np.radians(k * 360.0 / NSET)
        dy, dx = np.cos(a), np.sin(a)                 # verso sopravvento
        py, px = -dx, dy
        hr = raggio(dy, dx)
        dt = prima(hr)
        d15 = (r3 >= dt[None]) & (r3 <= dt[None] + 15.0)
        d3 = (r3 >= dt[None]) & (r3 <= dt[None] + 3.0)
        out["terra"][k] = dt
        out["h15"][k] = np.where(d15, hr, 0).max(axis=0)
        out["ang"][k] = np.degrees(np.arctan(hr / (r3 * 1000.0))).max(axis=0)
        out["hc"][k] = np.where(d3, hr, 0).sum(axis=0) / np.maximum(d3.sum(axis=0), 1)
        out["f5"][k] = (hr[passi <= 5] > 0).mean(axis=0)
        out["f15"][k] = (hr[passi <= 15] > 0).mean(axis=0)
        val = np.zeros(LA.shape); nval = np.zeros(LA.shape)
        for dd in (2.0, 4.0, 6.0):
            r = np.minimum(dt + dd, RMAX)
            cy = LA + r * dy / kmy; cx = LO + r * dx / kmx
            linea = quota(cy, cx)
            lati = np.max([quota(cy + q * py / kmy, cx + q * px / kmx) for q in (-4, -3, -2, 2, 3, 4)], axis=0)
            ok = dt < RMAX
            val += np.where(ok, lati - linea, 0); nval += ok
        out["valle"][k] = val / np.maximum(nval, 1)
        hs = raggio(-dy, -dx)                         # sottovento
        ds = prima(hs)
        d10 = (r3 >= ds[None]) & (r3 <= ds[None] + 10.0)
        out["sotto"][k] = ds
        out["hsotto"][k] = np.where(d10, hs, 0).max(axis=0)
        out["angsotto"][k] = np.degrees(np.arctan(hs / (r3 * 1000.0))).max(axis=0)
        # ai lati: le due perpendicolari (lsx, ldx); per il canale contano la somma e il minimo
        out["lsx"][k] = prima(raggio(py, px))
        out["ldx"][k] = prima(raggio(-py, -px))
        del hr, hs
    # ventaglio ±22,5°: i tre settori vicini
    out["angw"] = np.max([np.roll(out["ang"], s, axis=0) for s in (-1, 0, 1)], axis=0)
    out["fw15"] = np.mean([np.roll(out["f15"], s, axis=0) for s in (-1, 0, 1)], axis=0)
    # intorno, senza direzione
    rr = []
    for k in range(NSET):
        a = np.radians(k * 360.0 / NSET)
        rr.append(np.stack([quota(LA + r * np.cos(a) / kmy, LO + r * np.sin(a) / kmx) for r in passi[passi <= 10][::2]]))
    rr = np.concatenate(rr)
    fisse = {"h10": rr.max(axis=0), "hm10": np.where(rr > 0, rr, 0).sum(axis=0) / np.maximum((rr > 0).sum(axis=0), 1),
             "terra10": (rr > 0).mean(axis=0)}
    return out, fisse


def interp(f, da, j, i):
    """f: (16, ny, nx); da: gradi da cui viene il vento nelle celle (j, i) -> valore interpolato tra i due settori vicini."""
    x = da / (360.0 / NSET)
    k0 = np.floor(x).astype(int) % NSET
    w = (x - np.floor(x)).astype("f4")
    k1 = (k0 + 1) % NSET
    return (1 - w) * f[k0, j, i] + w * f[k1, j, i]


def main():
    cart, ril_p, usc = sys.argv[1:4]
    lens_cart = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] != "-" else None
    varianti = dict(a.split("=", 1) for a in sys.argv[5:])
    ril = Rilievo(ril_p)
    righe = {}
    nscena = 0

    def metti(k, v):
        righe.setdefault(k, []).append(np.asarray(v))

    for area in AREE:
        ia = AREE_TUTTE.index(area)
        z, raw = campi(cart, area)
        lat, lon = z["lat"], z["lon"]
        ns, _, ny, nx = raw.shape
        R = z["radar"].astype("f4")
        mu, mv = raw[:, 2], raw[:, 3]
        M = np.hypot(mu, mv)
        da = ((np.degrees(np.arctan2(-mu, -mv)) + 360) % 360).astype("f4")
        dist = z["dist"].astype("f4")
        D = np.broadcast_to(dist, R.shape)
        # la scena al largo
        largo = dist >= LARGO
        if largo.sum() < 20:
            largo = dist >= np.quantile(dist[dist > 0], 0.75)
        Ml = M[:, largo]; Rl = R[:, largo]
        S0 = np.nanmedian(Ml, axis=1)
        ul = np.nanmean(mu[:, largo], axis=1); vl = np.nanmean(mv[:, largo], axis=1)
        D0 = (np.degrees(np.arctan2(-ul, -vl)) + 360) % 360
        okl = np.isfinite(Rl) & np.isfinite(Ml)
        k0 = np.where(okl.sum(axis=1) >= 10, np.nansum(np.where(okl, Rl * Ml, 0), axis=1) /
                      np.maximum(np.nansum(np.where(okl, Ml * Ml, 0), axis=1), 1e-6), np.nan)
        nlargo = okl.sum(axis=1)
        st, fisse = statiche(lat, lon, ril)
        ok = np.isfinite(R) & np.isfinite(M) & (M > VMIN) & (D <= DMAX)
        si, sj, sk = np.nonzero(ok)
        print(f"{area}: {ns} passaggi, {len(si)} celle", flush=True)
        metti("area", np.full(len(si), ia, dtype="i1"))
        metti("scena", (nscena + si).astype("i4"))
        metti("j", sj.astype("i2")); metti("i", sk.astype("i2"))
        metti("radar", R[si, sj, sk]); metti("mod", M[si, sj, sk]); metti("dir", da[si, sj, sk])
        metti("dist", D[si, sj, sk])
        metti("cmems", z["cmems"][si, sj, sk].astype("f4"))
        Msh = np.hypot(raw[:, 0], raw[:, 1])
        metti("shift", Msh[si, sj, sk])
        for nome, f in st.items():
            metti(nome, interp(f, da[si, sj, sk], sj, sk).astype("f4"))
        for nome, f in fisse.items():
            metti(nome, f[sj, sk].astype("f4"))
        metti("S0", S0[si].astype("f4")); metti("D0", D0[si].astype("f4")); metti("k0", k0[si].astype("f4"))
        metti("nlargo", nlargo[si].astype("i4"))
        dd = (da[si, sj, sk] - D0[si] + 540) % 360 - 180
        metti("dDir", dd.astype("f4"))
        t = z["t"][si]
        metti("t", t.astype("i8"))
        metti("ora", ((t // 3600) % 24).astype("i1"))
        mese = (np.array(t, dtype="datetime64[s]").astype("datetime64[M]").astype(int) % 12 + 1).astype("i1")
        metti("mese", mese)
        if lens_cart:
            zl, rl_ = campi(lens_cart, area)
            assert np.array_equal(zl["t"], z["t"])
            L = np.hypot(rl_[:, 6], rl_[:, 7])
            metti("q_lens054", (L[si, sj, sk] / M[si, sj, sk]).astype("f4"))
        for nome, cv in varianti.items():
            if not os.path.exists(os.path.join(cv, f"sho_{area}.bin")):
                metti("q_" + nome, np.full(len(si), np.nan, dtype="f4"))
                continue
            zv, rv = campi(cv, area)
            assert np.array_equal(zv["t"], z["t"])
            L = np.hypot(rv[:, 6], rv[:, 7])
            metti("q_" + nome, (L[si, sj, sk] / M[si, sj, sk]).astype("f4"))
        nscena += ns
        del raw, R, M, da
    T = {k: np.concatenate(v) for k, v in righe.items()}
    np.savez(usc, **T)
    print(f"{len(T['radar'])} righe, {nscena} passaggi -> {usc}")


if __name__ == "__main__":
    main()
