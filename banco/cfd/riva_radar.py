"""Quanto costa, dove il radar vede, riaccendere il rilievo pieno vicino alla riva (LENS_RIVA "d1,d2").

  python3 riva_radar.py <tabella6.npz> <radar_vero_tabella6.npy> <cartella campi_<area>.npz> <uscita .json>
         area=atlante.lna [area=atlante.lna ...]

Per ogni area: la risposta dell'atlante di Lens (classe Fr 0,7) per 16 direzioni, media sulle celle di mare dentro ogni
cella del radar, divisa per la mediana sul mare a più di 6 km da terra. In ogni riga della tabella (cella e passaggio)
la direzione è quella del modello al largo. Lens 0.54 è la base (modello letto dal mare per q_lens054); la variante
moltiplica per (risposta dell'atlante) ^ w(distanza), con w = 1 entro d1 e 0 oltre d2. Nessuna stima: si guarda solo
quanto cambiano l'errore medio e l'errore sistematico per cella nelle fasce che il radar vede (1,5-2, 2-3 km).
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
from controlla_atlante import leggi   # noqa: E402

AREE14 = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova",
          "pontine", "corsica-nord", "monaco", "cicladi"]
VARIANTI = [(0.5, 2.0), (1.0, 2.5), (1.0, 3.0), (1.5, 3.0)]


def risposta_su_radar(p, rlat, rlon, fr=0.7):
    A = leggi(p)
    k = int(np.argmin(np.abs(np.array(A["fr"]) - fr)))
    ue, ve, un, vn = A["campi"][k]
    s, n, w, e = A["box"]
    lat = n - (np.arange(A["ny"]) + 0.5) * (n - s) / A["ny"]; lon = w + (np.arange(A["nx"]) + 0.5) * (e - w) / A["nx"]
    m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
    dist = distance_transform_edt(mare) * A["dx"] / 1000.0
    dla = abs(rlat[1] - rlat[0]); dlo = abs(rlon[1] - rlon[0])
    out = np.full((16, len(rlat), len(rlon)), np.nan)
    for kk, d in enumerate(np.arange(16) * 22.5):
        u0, v0 = -np.sin(np.radians(d)), -np.cos(np.radians(d))
        sp = np.hypot(u0 * ue + v0 * un, u0 * ve + v0 * vn)
        off = np.median(sp[(dist >= 6)]) if (dist >= 6).sum() > 50 else np.median(sp[mare])
        r = np.where(mare, sp / off, np.nan)
        for jr, a in enumerate(rlat):
            rj = np.where(np.abs(lat - a) <= dla / 2)[0]
            if not len(rj):
                continue
            for ir, b in enumerate(rlon):
                ri = np.where(np.abs(lon - b) <= dlo / 2)[0]
                if not len(ri):
                    continue
                blocco = r[rj[0]:rj[-1] + 1, ri[0]:ri[-1] + 1]
                if np.isfinite(blocco).mean() >= 0.3:
                    out[kk, jr, ir] = np.nanmean(blocco)
    return out


def peso(d, d1, d2):
    x = np.clip((d - d1) / (d2 - d1), 0, 1)
    return 0.5 + 0.5 * np.cos(np.pi * x)


def main():
    tab_p, vero_p, campi_d, usc = sys.argv[1:5]
    atl = dict(x.split("=", 1) for x in sys.argv[5:])
    T = np.load(tab_p, mmap_mode="r"); vero_all = np.load(vero_p, mmap_mode="r")
    area_all = np.asarray(T["area"])
    ris = {}
    for nome, p in atl.items():
        ia = AREE14.index(nome); m = area_all == ia
        C = np.load(os.path.join(campi_d, f"campi_{nome}.npz"))
        R = risposta_su_radar(p, C["lat"].astype(float), C["lon"].astype(float))
        j, i = T["j"][m].astype(int), T["i"][m].astype(int)
        vero = vero_all[m]; mod = T["mod"][m]; q = T["q_lens054"][m]; D0 = T["D0"][m]; dist = T["dist"][m]
        x = (D0 % 360) / 22.5; k0 = np.floor(x).astype(int) % 16; w = x - np.floor(x); k1 = (k0 + 1) % 16
        r = (1 - w) * R[k0, j, i] + w * R[k1, j, i]
        ok = np.isfinite(r) & np.isfinite(vero) & (mod > 2) & (q > 0) & (r > 0)
        base = mod * q
        cel = j * 1000 + i
        voce = {}
        for fascia, (lo, hi) in (("1.5-2", (1.5, 2.0)), ("2-3", (2.0, 3.0)), ("3-5", (3.0, 5.0))):
            f = ok & (dist >= lo) & (dist < hi)
            if f.sum() < 200:
                continue
            def misura(prev):
                e = (prev[f] - vero[f]) / 0.5144
                ids, inv = np.unique(cel[f], return_inverse=True); nn = np.bincount(inv); me = np.bincount(inv, e) / nn
                return dict(medio=round(float(e.mean()), 3), quadratico=round(float(np.sqrt(np.mean(e ** 2))), 3),
                            sistematico=round(float(np.sqrt(np.sum(nn * me ** 2) / nn.sum())), 3))
            v = {"righe": int(f.sum()), "lens054": misura(base), "rilievo_pieno": misura(base * r)}
            for d1, d2 in VARIANTI:
                v[f"riva_{d1:g}_{d2:g}"] = misura(base * r ** peso(dist, d1, d2))
            voce[fascia] = v
        ris[nome] = voce
        print(nome, json.dumps({f: {k: (vv["sistematico"], vv["medio"]) if isinstance(vv, dict) else vv for k, vv in v.items()} for f, v in voce.items()}), flush=True)
    json.dump(ris, open(usc, "w"), indent=1)


if __name__ == "__main__":
    main()
