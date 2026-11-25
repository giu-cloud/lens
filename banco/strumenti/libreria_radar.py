"""Lens 2: la libreria degli analoghi del radar per un'area e un modello (formato LRA1) e il calcolo di riferimento.

  VERO=<radar_vero.npy> python3 libreria_radar.py scrivi <tabella.npz> <cartella del metro> <modello> <uscita> [aree] [--cassaforte]
  VERO=<radar_vero.npy> python3 libreria_radar.py prova  <tabella.npz> <cartella librerie>   (giudica le librerie sulla cassaforte)
  VERO=<radar_vero.npy> python3 libreria_radar.py casi   <tabella.npz> <cartella librerie> <uscita.json> (casi per le prove in Kotlin)

Con --cassaforte la libreria si costruisce senza le settimane della cassaforte (per giudicarla); senza, con tutto.

Formato LRA1 (gzip, little-endian):
  "LRA1", int32 versione (1)
  area, modello: int16 lunghezza + UTF-8
  float64 lat0, lon0, dlat, dlon; int32 ny, nx           griglia del radar (centri delle celle)
  int32 nc; int32[nc] cella (j*nx + i); float32[nc] distanza dalla terra (km)
  float32 K, k_an, sd_v, sd_925, sd_ora, vmin, wbar        parametri
  float32[5*nc] armoniche (c, a1, b1, a2, b2; blocco per coefficiente)
  int32 ns; int64[ns] t; float32[ns] a0, b0, a925, b925; uint8[ns] sera
      a0 = S0·sin(D0), b0 = S0·cos(D0): vento al largo del modello (S0 forza mediana, D0 direzione di provenienza)
      a925, b925 lo stesso a 925 hPa; sera = 1 per i passaggi dopo le 12 UTC
  uint8[ns*nc] M: vento del modello nella cella, decimi di m/s (255 = manca)
  uint8[ns*nc] Q: rapporto vento vero (radar) ÷ modello, round(64·ln Q + 128) (255 = manca)
Il fattore per un'ora di previsione: le K scene più vicine nello spazio (a0, b0)/sd_v, (a925, b925)/sd_925, sera/sd_ora
(ordinate in modo stabile), pesi p = exp(−d²/2 mediana(d²)), poi per cella
  G = (Σ p·Q·M² + k_an·wbar·base) / (Σ p·M² + k_an·wbar),  base = armoniche sulla direzione del modello nella cella,
con solo le celle e scene dove M > vmin.
"""
import gzip
import io
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confronto_metodi import AREE, carica, misure
from metodi_area import PAR, Area, _armoniche

PARAMETRI = dict(PAR, sd_st=1000.0, sd_925=3.0)


def _s(b, x):
    e = x.encode("utf-8"); b.write(struct.pack("<h", len(e))); b.write(e)


def scrivi(path, area, modello, z, A, tr, p):
    coef = _armoniche(A, tr, p)
    R, M = A.R[tr], A.M[tr]
    ok = np.isfinite(R) & np.isfinite(M) & (M > p["vmin_an"])
    wbar = float(np.sum(np.where(ok, M * M, 0)) / max(ok.sum(), 1))
    celle = np.nonzero(np.isfinite(A.M).any(axis=0))[0]
    lat, lon = z["lat"], z["lon"]
    dist = z["dist"].ravel()[celle]
    Mq = np.where(np.isfinite(M[:, celle]), np.clip(np.round(M[:, celle] / 0.1), 0, 254), 255).astype(np.uint8)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.log(np.clip(R[:, celle] / M[:, celle], np.exp(-1.98), np.exp(1.98)))
    Qq = np.where(np.isfinite(q) & np.isfinite(M[:, celle]), np.clip(np.round(64 * q + 128), 0, 254), 255).astype(np.uint8)
    a = np.radians(A.D0[tr]); b = np.radians(A.d925[tr])
    b_ = io.BytesIO()
    b_.write(b"LRA1"); b_.write(struct.pack("<i", 1)); _s(b_, area); _s(b_, modello)
    b_.write(struct.pack("<4d2i", float(lat[0]), float(lon[0]), float(lat[1] - lat[0]), float(lon[1] - lon[0]), len(lat), len(lon)))
    b_.write(struct.pack("<i", len(celle))); b_.write(celle.astype("<i4").tobytes()); b_.write(dist.astype("<f4").tobytes())
    b_.write(struct.pack("<7f", p["K"], p["k_an"], p["sd_v"], p["sd_925"], p["sd_ora"], p["vmin_an"], wbar))
    b_.write(coef[:, celle].astype("<f4").tobytes())
    b_.write(struct.pack("<i", len(tr))); b_.write(A.t[tr].astype("<i8").tobytes())
    for arr in (A.S0[tr] * np.sin(a), A.S0[tr] * np.cos(a), A.v925[tr] * np.sin(b), A.v925[tr] * np.cos(b)):
        b_.write(arr.astype("<f4").tobytes())
    b_.write(A.sera[tr].astype(np.uint8).tobytes())
    b_.write(Mq.tobytes()); b_.write(Qq.tobytes())
    with gzip.open(path, "wb", compresslevel=9) as f:
        f.write(b_.getvalue())
    return len(celle), len(tr)


class Libreria:
    def __init__(self, path):
        d = gzip.open(path, "rb").read(); o = 0

        def rd(fmt):
            nonlocal o
            v = struct.unpack_from(fmt, d, o); o += struct.calcsize(fmt); return v

        def rs():
            n, = rd("<h"); nonlocal o
            s = d[o:o + n].decode("utf-8"); o += n; return s

        def arr(dt, n):
            nonlocal o
            a = np.frombuffer(d, dtype=dt, count=n, offset=o); o += a.nbytes; return a
        assert d[:4] == b"LRA1"; o = 4
        self.versione, = rd("<i"); self.area = rs(); self.modello = rs()
        self.lat0, self.lon0, self.dlat, self.dlon, self.ny, self.nx = rd("<4d2i")
        self.nc, = rd("<i"); self.celle = arr("<i4", self.nc); self.dist = arr("<f4", self.nc)
        self.K, self.k_an, self.sd_v, self.sd_925, self.sd_ora, self.vmin, self.wbar = rd("<7f")
        self.coef = arr("<f4", 5 * self.nc).reshape(5, self.nc)
        self.ns, = rd("<i"); self.t = arr("<i8", self.ns)
        self.a0, self.b0, self.a925, self.b925 = (arr("<f4", self.ns) for _ in range(4))
        self.sera = arr("u1", self.ns)
        Mq = arr("u1", self.ns * self.nc).reshape(self.ns, self.nc); Qq = arr("u1", self.ns * self.nc).reshape(self.ns, self.nc)
        self.M = np.where(Mq < 255, Mq.astype("f4") * 0.1, np.nan)
        self.Q = np.where(Qq < 255, np.exp((Qq.astype("f4") - 128) / 64.0), np.nan)
        ok = np.isfinite(self.M) & np.isfinite(self.Q) & (self.M > self.vmin)
        self.ww = np.where(ok, self.M * self.M, 0.0).astype("f8")
        self.wy = np.where(ok, self.Q * self.M * self.M, 0.0).astype("f8")
        self.X = np.stack([self.a0 / self.sd_v, self.b0 / self.sd_v, self.sera / self.sd_ora,
                           self.a925 / self.sd_925, self.b925 / self.sd_925], axis=1).astype("f8")

    def fattore(self, a0, b0, a925, b925, sera, dir_celle, dispersione=False):
        """Fattore per cella (nc) per un'ora: dati della situazione e direzione del modello nelle celle (gradi).
        Con dispersione=True anche la dispersione degli analoghi attorno al fattore (scarto pesato del rapporto):
        misura quanto fidarsi, e l'errore relativo che resta le va dietro (prova del 25/11)."""
        x = np.array([a0 / self.sd_v, b0 / self.sd_v, sera / self.sd_ora, a925 / self.sd_925, b925 / self.sd_925])
        d2 = np.sum((self.X - x) ** 2, axis=1)
        K = min(int(self.K), self.ns)
        idx = np.argsort(d2, kind="stable")[:K]
        med = np.median(d2[idx])
        pk = np.exp(-0.5 * d2[idx] / max(med, 1e-6))
        th = np.radians(dir_celle)
        base = (self.coef[0] + self.coef[1] * np.cos(th) + self.coef[2] * np.sin(th) + self.coef[3] * np.cos(2 * th)
                + self.coef[4] * np.sin(2 * th))
        lam = self.k_an * self.wbar
        G = (pk @ self.wy[idx] + lam * base) / (pk @ self.ww[idx] + lam)
        if not dispersione:
            return G
        q = np.where(self.ww[idx] > 0, self.wy[idx] / np.where(self.ww[idx] > 0, self.ww[idx], 1), 0.0)
        w = pk[:, None] * self.ww[idx]
        sp = np.sqrt((w * (q - G) ** 2).sum(0) / np.maximum(w.sum(0), 1e-9))
        return G, np.where(w.sum(0) > 0, sp, np.nan)


def _aree(T, nomi):
    return [AREE.index(a) for a in nomi] if nomi else sorted(set(T["area"].tolist()))


def main():
    cmd = sys.argv[1]
    if cmd == "scrivi":
        tab, metro, modello, usc = sys.argv[2:6]
        resto = sys.argv[6:]
        cassa = "--cassaforte" in resto
        nomi = [a for a in resto if a in AREE]
        T = carica(tab); os.makedirs(usc, exist_ok=True)
        for ia in _aree(T, nomi):
            area = AREE[ia]
            A = Area(T, ia, np.nonzero(T["area"] == ia)[0])
            z = np.load(os.path.join(metro, f"campi_{area}.npz"))
            A.ny, A.nx = len(z["lat"]), len(z["lon"])     # la griglia intera del metro
            A2 = Area.__new__(Area); A2.__dict__.update(A.__dict__)
            # ricostruisce le matrici sulla griglia intera
            nc = A.ny * A.nx; righe = A.r
            cella = T["j"][righe].astype(int) * A.nx + T["i"][righe].astype(int)
            for nome, col in (("R", "radar"), ("M", "mod"), ("D", "dir")):
                X = np.full((len(A.scene), nc), np.nan, "f4"); X[A.inv, cella] = T[col][righe]; setattr(A2, nome, X)
            A2.cella = cella
            tr = np.nonzero(~A2.cassa)[0] if cassa else np.arange(len(A2.scene))
            p = os.path.join(usc, f"lra1_{modello}_{area}.lra")
            nc_, ns_ = scrivi(p, area, modello, z, A2, tr, PARAMETRI)
            print(f"{area}: {nc_} celle, {ns_} passaggi, {os.path.getsize(p) / 1e6:.2f} MB", flush=True)
    elif cmd in ("prova", "casi"):
        tab, cart = sys.argv[2:4]
        T = carica(tab)
        pred = np.full(len(T["radar"]), np.nan, "f4")
        casi = []
        for f in sorted(os.listdir(cart)):
            if not f.endswith(".lra"):
                continue
            L = Libreria(os.path.join(cart, f))
            ia = AREE.index(L.area)
            righe = np.nonzero((T["area"] == ia) & T["cassa"])[0]
            cel = T["j"][righe].astype(int) * L.nx + T["i"][righe].astype(int)
            pos = {c: k for k, c in enumerate(L.celle)}
            for s in np.unique(T["scena"][righe]):
                rr = righe[T["scena"][righe] == s]; cc = cel[T["scena"][righe] == s]
                r0 = rr[0]
                a = np.radians(T["D0"][r0]); b = np.radians(np.nan_to_num(T["d925"][r0]))
                f_ = (float(T["S0"][r0] * np.sin(a)), float(T["S0"][r0] * np.cos(a)),
                      float(np.nan_to_num(T["v925"][r0]) * np.sin(b)), float(np.nan_to_num(T["v925"][r0]) * np.cos(b)),
                      float(T["sera"][r0]))
                dirc = np.zeros(L.nc, "f4")
                k = np.array([pos.get(c, -1) for c in cc]); ok = k >= 0
                dirc[k[ok]] = T["dir"][rr[ok]]
                G, SP = L.fattore(*f_, dirc, dispersione=True)
                pred[rr[ok]] = G[k[ok]]
                if cmd == "casi" and len([c for c in casi if c["area"] == L.area]) < 3:
                    sel = np.nonzero(ok)[0][:: max(1, ok.sum() // 40)][:40]
                    casi.append({"libreria": f, "area": L.area, "dati": f_, "celle": [int(k[i]) for i in sel],
                                 "direzioni": [float(dirc[k[i]]) for i in sel], "atteso": [float(G[k[i]]) for i in sel],
                                 "dispersione": [float(SP[k[i]]) for i in sel],
                                 "dir_tutte": [round(float(x), 2) for x in dirc]})
        ix = np.nonzero(np.isfinite(pred))[0]
        r = misure(T, ix, pred[ix], "ass"); r0 = misure(T, ix, np.ones(len(ix), "f4"), "ass")
        print(f"modello   rmse {r0['rmse']:.3f} sistematico {r0['sistematico']:.3f} (0-3 km {r0['sis0-3']:.3f})")
        print(f"libreria  rmse {r['rmse']:.3f} sistematico {r['sistematico']:.3f} (0-3 km {r['sis0-3']:.3f})")
        if cmd == "casi":
            json.dump(casi, open(sys.argv[4], "w"))
            print(f"{len(casi)} casi -> {sys.argv[4]}")


if __name__ == "__main__":
    main()
