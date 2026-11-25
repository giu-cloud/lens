"""La CFD contro le stazioni sulla costa di un sito: forza e direzione, settore per settore, a tutte le ore.

  python3 stazioni_cfd.py --sito stretto --cfd <cartella cfd_<sito>_*.npz> --mare <terreno .npz o atlante .lna>
        --modello <modello_siti.npz> --elenco <stazioni_elenco.json> --stazioni <cartella stazioni> --uscita <.json>
        [--atlante <atlante Lens .lna dell'area>] [--fr 1.0] [--vmin 5]

Il radar non vede il primo chilometro e mezzo dalla riva e non misura la direzione: qui giudicano le stazioni.
Per ogni stazione dentro il dominio della CFD e per ogni ora (ogni 3 ore, ECMWF al largo del sito: S0, D0):
  - forza: rapporto stazione / S0 contro il rapporto della CFD nel punto della stazione (vento a 10 m sul suolo,
    diviso la mediana della CFD sul mare a più di 6 km da terra), interpolato tra le due direzioni vicine a D0;
  - direzione: la stazione meno D0 contro la rotazione della CFD nel punto.
Si usano solo le ore con S0 >= vmin (5 m/s: meno brezze) e la stazione sopra 1 m/s.
La stazione ha un'esposizione sua (altezza, alberi, edifici) che nessun modello a 150 m conosce: per la forza si
giudica il disegno per settori (correlazione dei logaritmi tra i settori) e la pendenza ora per ora con una costante
per stazione stimata fuori dalla cassaforte; per la direzione l'errore medio con e senza la rotazione della CFD.
Se c'è l'atlante di Lens dell'area, lo stesso per la risposta dell'atlante (conservazione della massa, classe --fr).
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
from taratura_radar import pila_cfd, mare_su   # noqa: E402

SETTIMANA = 7 * 86400


def bilineare(F, lat, lon, la, lo):
    """F (..., ny, nx) sulla griglia (lat da nord, lon); NaN se un vertice è NaN."""
    fy = (lat[0] - la) / (lat[0] - lat[1]); fx = (lo - lon[0]) / (lon[1] - lon[0])
    j0 = int(np.floor(fy)); i0 = int(np.floor(fx)); ty = fy - j0; tx = fx - i0
    if j0 < 0 or i0 < 0 or j0 + 1 >= len(lat) or i0 + 1 >= len(lon):
        return np.full(F.shape[:-2], np.nan)
    return ((1 - ty) * (1 - tx) * F[..., j0, i0] + (1 - ty) * tx * F[..., j0, i0 + 1] + ty * (1 - tx) * F[..., j0 + 1, i0]
            + ty * tx * F[..., j0 + 1, i0 + 1])


def vicini(da, theta):
    """Indici delle due direzioni calcolate attorno a theta e peso della seconda (anche con direzioni mancanti)."""
    da = np.asarray(da, float); th = np.asarray(theta, float) % 360.0
    k1 = np.searchsorted(da, th, side="right") % len(da)
    k0 = (k1 - 1) % len(da)
    gap = (da[k1] - da[k0]) % 360.0
    gap = np.where(gap == 0, 360.0, gap)
    w = ((th - da[k0]) % 360.0) / gap
    return k0, k1, w


def interp_dir(da, Y, theta):
    """Y (k, n) o (k,) alle direzioni theta (n,): lineare tra le due direzioni calcolate vicine."""
    k0, k1, w = vicini(da, theta)
    return (1 - w) * Y[k0] + w * Y[k1]


def ang(a):
    return (np.asarray(a) + 540.0) % 360.0 - 180.0


def risposta_atlante(p, fr, la, lo):
    """Rapporto e rotazione della risposta dell'atlante di Lens nel punto, per le 16 direzioni (classe fr)."""
    from controlla_atlante import leggi
    A = leggi(p)
    k = int(np.argmin(np.abs(np.array(A["fr"]) - fr)))
    ue, ve, un, vn = A["campi"][k]
    s, n, w, e = A["box"]
    lat = n - (np.arange(A["ny"]) + 0.5) * (n - s) / A["ny"]; lon = w + (np.arange(A["nx"]) + 0.5) * (e - w) / A["nx"]
    m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
    dist = distance_transform_edt(mare) * A["dx"] / 1000.0
    da = np.arange(16) * 22.5
    rap, rot = [], []
    for d in da:
        u0, v0 = -np.sin(np.radians(d)), -np.cos(np.radians(d))
        U = u0 * ue + v0 * un; V = u0 * ve + v0 * vn
        sp = np.hypot(U, V)
        off = np.median(sp[dist >= 6]) if (dist >= 6).sum() > 20 else np.median(sp[mare])
        uu, vv = bilineare(np.stack([U, V]), lat, lon, la, lo)
        rap.append(np.hypot(uu, vv) / off)
        rot.append(ang((np.degrees(np.arctan2(-uu, -vv)) + 360) % 360 - d))
    return da, np.array(rap), np.array(rot)


def main():
    ap = argparse.ArgumentParser()
    for k in ("sito", "cfd", "mare", "modello", "elenco", "stazioni", "uscita"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--atlante", default=""); ap.add_argument("--fr", type=float, default=1.0)
    ap.add_argument("--vmin", type=float, default=5.0)
    a = ap.parse_args()
    da, U, V, lat, lon = pila_cfd(a.cfd, a.sito)
    S = np.hypot(U, V)
    mare = mare_su(a.mare, lat, lon)
    dx_km = abs(lat[1] - lat[0]) * 111.19
    dist = distance_transform_edt(mare) * dx_km
    Soff = np.array([np.nanmedian(np.where((dist >= 6) & mare, S[k], np.nan)) for k in range(len(da))])
    M = np.load(a.modello)
    tm, S0, D0 = M[f"{a.sito}_t"], M[f"{a.sito}_S0"], M[f"{a.sito}_D0"]
    la0, lo0 = float(np.mean(lat)), float(np.mean(lon))
    elenco = json.load(open(a.elenco))
    ris = {"sito": a.sito, "direzioni": da.tolist(), "vmin": a.vmin, "stazioni": []}
    for st in elenco:
        p = os.path.join(a.stazioni, f"{st['fonte']}_{st['codice']}.npz")
        if not os.path.exists(p):
            continue
        uv = bilineare(np.stack([U, V]), lat, lon, st["lat"], st["lon"])     # (2, k)
        if not np.isfinite(uv).all():
            continue
        rc = np.hypot(uv[0], uv[1]) / Soff
        rotc = ang((np.degrees(np.arctan2(-uv[0], -uv[1])) + 360) % 360 - da)
        Z = np.load(p); o = np.argsort(Z["t"]); t, ve, di = Z["t"][o], Z["vento"][o].astype(float), Z["direzione"][o].astype(float)
        k = np.searchsorted(t, tm); k = np.clip(k, 1, len(t) - 1)
        kk = np.where(np.abs(t[k - 1] - tm) < np.abs(t[k] - tm), k - 1, k)
        vicino = np.abs(t[kk] - tm) <= 20 * 60
        oss_v, oss_d = ve[kk], di[kk]
        ok = vicino & (S0 >= a.vmin) & np.isfinite(oss_v) & (oss_v > 1.0) & np.isfinite(oss_d)
        if ok.sum() < 100:
            continue
        cass = (tm // SETTIMANA) % 5 == 2
        th = D0[ok]
        y = np.log(oss_v[ok] / S0[ok])
        pc = np.log(interp_dir(da, rc, th)); rotp = interp_dir(da, rotc, th)
        dev = ang(oss_d[ok] - th)
        voce = dict(nome=st["nome"], fonte=st["fonte"], codice=st["codice"], lat=st["lat"], lon=st["lon"], quota=st.get("quota"),
                    km_dal_centro=round(float(np.hypot((st["lat"] - la0) * 111.19, (st["lon"] - lo0) * 111.19 * np.cos(np.radians(la0)))), 1),
                    ore=int(ok.sum()), cfd_rapporto=np.round(rc, 3).tolist(), cfd_rotazione=np.round(rotc, 1).tolist())
        modelli = {"cfd": (pc, rotp)}
        if a.atlante:
            dA, rA, rotA = risposta_atlante(a.atlante, a.fr, st["lat"], st["lon"])
            if np.isfinite(rA).all():
                modelli["atlante"] = (np.log(interp_dir(dA, rA, th)), interp_dir(dA, rotA, th))
                voce["atlante_rapporto"] = np.round(rA, 3).tolist(); voce["atlante_rotazione"] = np.round(rotA, 1).tolist()
        # settori: mediana osservata contro mediana prevista
        sett = (np.round(th / 22.5).astype(int)) % 16
        tab = []
        for s_ in range(16):
            m = sett == s_
            if m.sum() < 15:
                continue
            riga = dict(settore=s_ * 22.5, ore=int(m.sum()), oss_log=float(np.median(y[m])),
                        oss_rot=float(np.degrees(np.arctan2(np.mean(np.sin(np.radians(dev[m]))), np.mean(np.cos(np.radians(dev[m])))))))
            for nome, (pl, pr) in modelli.items():
                riga[nome + "_log"] = float(np.median(pl[m])); riga[nome + "_rot"] = float(np.median(pr[m]))
            tab.append(riga)
        voce["settori"] = tab
        w = np.array([r["ore"] for r in tab], float)
        def cpes(x, z):
            x, z = np.asarray(x), np.asarray(z)
            if len(x) < 4 or np.std(x) == 0 or np.std(z) == 0:
                return None
            xm, zm = np.average(x, weights=w), np.average(z, weights=w)
            return float(np.sum(w * (x - xm) * (z - zm)) / np.sqrt(np.sum(w * (x - xm) ** 2) * np.sum(w * (z - zm) ** 2)))
        voce["giudizio"] = {}
        tr, te = ~cass[ok], cass[ok]
        oss_rot = [r["oss_rot"] for r in tab]
        voce["giudizio"]["direzione_err_modello"] = float(np.mean(np.abs(dev[te])))
        for nome, (pl, pr) in modelli.items():
            g = {"disegno_settori": cpes([r["oss_log"] for r in tab], [r[nome + "_log"] for r in tab]),
                 "rotazione_settori": cpes(oss_rot, [r[nome + "_rot"] for r in tab])}
            # ora per ora: y = c + b * previsto (stimati fuori dalla cassaforte), errore nella cassaforte
            X = np.column_stack([np.ones(tr.sum()), pl[tr]]); c = np.linalg.lstsq(X, y[tr], rcond=None)[0]
            e_mod = y[te] - y[tr].mean(); e_p = y[te] - (c[0] + c[1] * pl[te])
            g.update(pendenza=float(c[1]), rmse_log_costante=float(np.sqrt(np.mean(e_mod ** 2))), rmse_log_previsto=float(np.sqrt(np.mean(e_p ** 2))),
                     direzione_err_con_rotazione=float(np.mean(np.abs(ang(dev[te] - pr[te])))),
                     direzione_err_per_peso={str(k): float(np.mean(np.abs(ang(dev[te] - k * pr[te])))) for k in (0.5, 1.0, 1.5, 2.0)},
                     direzione_err_vento_forte=float(np.mean(np.abs(ang(dev[te & (oss_v[ok] > 5)] - pr[te & (oss_v[ok] > 5)])))) if (te & (oss_v[ok] > 5)).sum() > 30 else None,
                     direzione_err_vento_forte_modello=float(np.mean(np.abs(dev[te & (oss_v[ok] > 5)]))) if (te & (oss_v[ok] > 5)).sum() > 30 else None)
            voce["giudizio"][nome] = g
        ris["stazioni"].append(voce)
        gg = voce["giudizio"]
        txt = " · ".join(f"{n}: disegno {gg[n]['disegno_settori'] if gg[n]['disegno_settori'] is None else round(gg[n]['disegno_settori'], 2)}, "
                         f"pendenza {gg[n]['pendenza']:.2f}, log {gg[n]['rmse_log_costante']:.3f}->{gg[n]['rmse_log_previsto']:.3f}, "
                         f"dir {gg['direzione_err_modello']:.0f}->{gg[n]['direzione_err_con_rotazione']:.0f}°" for n in modelli)
        print(f"{st['nome'][:22]:22s} ({st['fonte']}, {st.get('quota')} m, {voce['ore']} ore) {txt}")
    json.dump(ris, open(a.uscita, "w"), indent=1)


if __name__ == "__main__":
    main()
