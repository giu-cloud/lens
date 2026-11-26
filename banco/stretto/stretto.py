"""Studio dello Stretto di Messina (26/11/2026): che cosa vede il radar e che cosa fanno i metodi.

Regimi dal vento del modello al largo (D0): N (315-45), S (135-225), altri. Per ogni regime:
  - carte composite del rapporto vento vero (radar) / modello, cella per cella (minimi quadrati: Σ vero·mod / Σ mod²);
  - le stesse per le previsioni: 0.54, A (primo Lens), C (A con il livello del radar a 30 km), B (fattore della
    libreria cella per cella), CFD (risposta neutra a 150 m, interpolata sulla direzione);
  - profilo lungo l'asse dello Stretto;
  - errori sulle settimane della cassaforte (libreria costruita senza di esse) per regime e distanza dalla terra.

  python3 stretto.py <uscita .npz> <uscita .json>
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates

qui = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(qui, "..", "strumenti"), os.path.join(qui, "..", "cfd")]
from controlla_atlante import leggi as leggi_lna      # noqa: E402
from libreria_radar import Libreria                    # noqa: E402
from stazioni_cfd import vicini                        # noqa: E402
from taratura_radar import pila_cfd, su_celle_radar    # noqa: E402
from lens2_su_primo import risposta, liscia            # noqa: E402

L_ = "/home/claude/lavoro"
ATL = "/home/claude/dati/shift/atlante_messina.lna"
CAMPI = f"{L_}/v_scia/campi_messina.npz"
LIBC = f"{L_}/libc6/lra1_ecmwf_messina.lra"
LIB = f"{L_}/librerie/lra1_ecmwf_messina.lra"
P07, P20 = f"{L_}/primo16/messina_fr0.7.bin", f"{L_}/primo16/messina_fr2.0.bin"
CFD = "/home/claude/dati/cfd/stretto"
SETTIMANA = 7 * 86400
DA = np.arange(16) * 22.5
REGIMI = {"N": (315, 45), "S": (135, 225), "altri": None}
# l'asse dello Stretto, da nord (Tirreno) a sud (Ionio)
ASSE = [(38.33, 15.70), (38.29, 15.67), (38.262, 15.645), (38.245, 15.63), (38.225, 15.62), (38.20, 15.615),
        (38.175, 15.605), (38.15, 15.60), (38.12, 15.59), (38.09, 15.585), (38.06, 15.575), (38.03, 15.565), (38.00, 15.56), (37.97, 15.555)]


def regime(d0):
    out = np.full(len(d0), "altri", dtype=object)
    out[(d0 >= 315) | (d0 < 45)] = "N"
    out[(d0 >= 135) & (d0 < 225)] = "S"
    return out


def main():
    usc_npz, usc_json = sys.argv[1:3]
    A = leggi_lna(ATL); m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
    C = np.load(CAMPI); rlat, rlon, dist_c = C["lat"].astype(float), C["lon"].astype(float), C["dist"]
    dkm = abs(rlat[1] - rlat[0]) * 111.19
    R07 = risposta(P07, mare, rlat, rlon, dist_c); R20 = risposta(P20, mare, rlat, rlon, dist_c)
    # 0.54 (il primo campo dei file di Campi, con la direzione del rilievo: la forza è quella della 0.54)
    from leggi_campi import leggi
    la, lo, F = leggi(P20, 16)
    S54 = np.stack([np.hypot(*F[k][0]) for k in range(16)])
    S54r, _ = su_celle_radar(S54, mare, la, lo, rlat, rlon)
    off = np.array([np.nanmedian(S54r[k][dist_c > 6]) for k in range(16)])
    R54 = S54r / off[:, None, None]
    # CFD: risposta neutra sulle celle del radar
    dac, Uc, Vc, clat, clon = pila_cfd(CFD, "stretto")
    Sc = np.hypot(Uc, Vc)
    mc = np.isfinite(Sc).all(0)
    Scr, _ = su_celle_radar(np.where(mc, Sc, 0.0), mc, clat, clon, rlat, rlon)
    offc = np.array([np.nanmedian(Scr[k][dist_c > 6]) for k in range(len(dac))])
    RC = Scr / offc[:, None, None]
    valido = np.isfinite(R07).all(0) & np.isfinite(R20).all(0)
    sig = 30.0 / 2 / dkm
    R07s = np.stack([liscia(R07[k], valido, sig) for k in range(16)]); R20s = np.stack([liscia(R20[k], valido, sig) for k in range(16)])

    T = np.load(f"{L_}/tabella6.npz", mmap_mode="r")
    sel = np.asarray(T["area"]) == 3
    g = lambda k: np.asarray(T[k][sel])                                     # noqa: E731
    vero = np.load(f"{L_}/radar_vero_tabella6.npy", mmap_mode="r")[sel]
    j, i = g("j").astype(int), g("i").astype(int)
    mod, D0, S0, dist, t, sc = g("mod"), g("D0"), g("S0"), g("dist"), g("t"), g("scena")
    qL, dT925, dirc, v925, d925, ora = g("q_lens054"), g("dT925"), g("dir"), g("v925"), g("d925"), g("ora")
    cass = ((t // SETTIMANA) % 5 == 2)
    neu = np.clip((4 - dT925) / 4, 0, 1); neu = np.where(np.isfinite(neu), neu, 0.5)
    k0, k1, w = vicini(DA, D0)
    su = lambda R: (1 - w) * R[k0, j, i] + w * R[k1, j, i]                  # noqa: E731
    rA = (1 - neu) * su(R07) + neu * su(R20)
    rAs = (1 - neu) * su(R07s) + neu * su(R20s)
    r54 = su(R54)
    kc0, kc1, wc = vicini(dac, D0)
    rC = (1 - wc) * RC[kc0, j, i] + wc * RC[kc1, j, i]

    # Lens 2: fattore della libreria senza cassaforte sulle scene della cassaforte; liscio a 30 km per C
    def fattori(p_lib, scene):
        Lb = Libreria(p_lib); pos = {c: k for k, c in enumerate(Lb.celle)}
        jl, il = np.array(Lb.celle) // Lb.nx, np.array(Lb.celle) % Lb.nx
        mk = np.zeros((Lb.ny, Lb.nx), bool); mk[jl, il] = True
        G = np.full(len(j), np.nan); G30 = np.full(len(j), np.nan)
        cella = j * Lb.nx + i
        for s_id in scene:
            r = np.where(sc == s_id)[0]
            a = np.radians(D0[r[0]]); b = np.radians(d925[r[0]])
            if not (np.isfinite(a) and np.isfinite(b) and np.isfinite(v925[r[0]])):
                continue
            dc = np.full(Lb.nc, D0[r[0]]); idx = np.array([pos.get(c, -1) for c in cella[r]])
            dc[idx[idx >= 0]] = dirc[r[idx >= 0]]
            gg = Lb.fattore(S0[r[0]] * np.sin(a), S0[r[0]] * np.cos(a), v925[r[0]] * np.sin(b), v925[r[0]] * np.cos(b),
                            1.0 if ora[r[0]] >= 12 else 0.0, dc)
            G[r[idx >= 0]] = gg[idx[idx >= 0]]
            gr = np.full((Lb.ny, Lb.nx), np.nan); gr[jl, il] = np.log(gg)
            ls = liscia(gr, mk & np.isfinite(gr), sig)
            G30[r] = np.exp(ls[j[r], i[r]])
        return G, G30
    Gc, G30c = fattori(LIBC, np.unique(sc[cass]))

    reg = regime(D0)
    ok = np.isfinite(vero) & (vero > 0.5) & (mod > 2.0) & np.isfinite(rA) & np.isfinite(rC) & (S0 > 2.0)
    prev = {"modello": mod, "0.54": mod * np.where(np.isfinite(qL) & (qL > 0), qL, r54), "A": mod * rA,
            "C": mod * rA * G30c / rAs, "B": mod * Gc, "CFD": mod * rC}
    cel = j * 100 + i
    ris = {"passaggi": {}, "errori": {}}

    def sist(e, mk):
        ids, inv = np.unique(cel[mk], return_inverse=True); n = np.bincount(inv); me = np.bincount(inv, e[mk]) / n
        return float(np.sqrt(np.sum(n * me ** 2) / n.sum()))
    for rg in ("N", "S", "altri", "tutti"):
        base = ok & cass & np.isfinite(Gc) & np.isfinite(G30c) & ((reg == rg) if rg != "tutti" else True)
        ris["passaggi"][rg] = int(len(np.unique(sc[base])))
        for fascia, f in (("tutto", base), ("entro 3 km", base & (dist < 3)), ("3-10 km", base & (dist >= 3) & (dist < 10))):
            if f.sum() < 100:
                continue
            for nome, p in prev.items():
                e = (p - vero) / 0.5144
                ris["errori"].setdefault(rg, {}).setdefault(fascia, {})[nome] = dict(
                    rmse=round(float(np.sqrt(np.mean(e[f] ** 2))), 2), media=round(float(np.mean(e[f])), 2), sist=round(sist(e, f), 2), righe=int(f.sum()))
    # carte composite su tutte le scene (radar, modello, A, CFD, 0.54) e sulla cassaforte (B, C)
    ny, nx = len(rlat), len(rlon)
    carte = {}
    for rg in ("N", "S"):
        for nome, p, msk in (("radar", vero, ok), ("A", mod * rA, ok), ("CFD", mod * rC, ok), ("0.54", prev["0.54"], ok),
                             ("B", prev["B"], ok & cass & np.isfinite(Gc)), ("C", prev["C"], ok & cass & np.isfinite(G30c)),
                             ("radar_cass", vero, ok & cass & np.isfinite(Gc))):
            f = msk & (reg == rg) & (S0 > 4)
            num = np.zeros(ny * nx); den = np.zeros(ny * nx); n = np.zeros(ny * nx)
            k = j[f] * nx + i[f]
            np.add.at(num, k, p[f] * mod[f]); np.add.at(den, k, mod[f] ** 2); np.add.at(n, k, 1)
            carte[f"{rg}_{nome}"] = np.where(n >= 5, num / np.maximum(den, 1e-9), np.nan).reshape(ny, nx)
            carte[f"{rg}_{nome}_n"] = n.reshape(ny, nx)
    # profilo lungo l'asse (interpolazione bilineare sulle carte, ogni 1 km)
    pts = []
    for (a1, b1), (a2, b2) in zip(ASSE[:-1], ASSE[1:]):
        dd = np.hypot((a2 - a1) * 111.19, (b2 - b1) * 111.19 * np.cos(np.radians(a1)))
        for s in np.arange(0, dd, 1.0):
            pts.append((a1 + (a2 - a1) * s / dd, b1 + (b2 - b1) * s / dd))
    pts.append(ASSE[-1]); pts = np.array(pts)
    fj = (pts[:, 0] - rlat[0]) / (rlat[1] - rlat[0]); fi = (pts[:, 1] - rlon[0]) / (rlon[1] - rlon[0])
    prof = {"km": list(np.round(np.r_[0, np.cumsum(np.hypot(np.diff(pts[:, 0]) * 111.19, np.diff(pts[:, 1]) * 111.19 * np.cos(np.radians(38.1))))], 2)),
            "lat": list(np.round(pts[:, 0], 4)), "lon": list(np.round(pts[:, 1], 4))}
    for kk, Z in carte.items():
        if kk.endswith("_n"):
            continue
        z = map_coordinates(np.nan_to_num(Z, nan=-9), [fj, fi], order=0, mode="nearest")
        prof[kk] = [None if v < -1 else round(float(v), 3) for v in z]
    ris["profilo"] = prof
    np.savez(usc_npz, lat=rlat, lon=rlon, dist=dist_c, **carte)
    json.dump(ris, open(usc_json, "w"), indent=1)
    for rg, d in ris["errori"].items():
        print(f"== {rg} ({ris['passaggi'][rg]} passaggi della cassaforte)")
        for fascia, dd in d.items():
            print(f"  {fascia:11s} " + "  ".join(f"{n}: {v['rmse']:.2f}/{v['media']:+.2f}/{v['sist']:.2f}" for n, v in dd.items()) + f"   ({next(iter(dd.values()))['righe']} righe)")


if __name__ == "__main__":
    main()
