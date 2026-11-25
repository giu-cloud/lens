"""Il vento di ECMWF al largo di ogni sito CFD, ogni 3 ore, per la prova con le stazioni.

  python3 modello_siti.py <cartella ecmwf passi 3-6> <cartella e912 passi 9-12> <uscita .npz> <rilievo di Shift .rgr>
         sito:LAT:LON[,sito:LAT:LON...]

Per ogni sito: i nodi di ECMWF entro 45 km dal centro (RAGGIO_KM) che stanno ad almeno 8 km dalla terra (sul rilievo di Shift); S0 = mediana della forza, D0 = direzione del vento medio (da dove viene).
Corse delle 00 e delle 12: passi 3, 6 (prima cartella) e 9, 12 (seconda), cioè tutte le 3 ore.
"""
import datetime as dt
import glob
import os
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strumenti"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stabilita import leggi          # noqa: E402
from rilievo import Rilievo          # noqa: E402


RAGGIO = float(os.environ.get("RAGGIO_KM", "45"))


def nodi_largo(lat, lon, ril, la0, lo0):
    """Nodi entro RAGGIO km dal centro senza terra entro 8 km (poi 6, 4) sul rilievo di Shift."""
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    kx = 111.19 * np.cos(np.radians(la0))
    vic = np.hypot((LA - la0) * 111.19, (LO - lo0) * kx) < RAGGIO
    dist = np.full(LA.shape, 0.0)
    for j, i in zip(*np.where(vic)):
        a, b = LA[j, i], LO[j, i]
        gla, glo, h = ril.griglia(0, a - 0.1, a + 0.1, b - 0.13, b + 0.13)
        GA, GO = np.meshgrid(gla, glo, indexing="ij")
        d = np.hypot((GA - a) * 111.19, (GO - b) * kx)
        dist[j, i] = d[h > 0].min() if (h > 0).any() else 11.0
    for soglia in (8, 6, 4):
        sel = vic & (dist >= soglia)
        if sel.sum() >= 3:
            return sel, soglia
    return vic & (dist >= 2), 2


def main():
    e36, e912, usc, ril_p = sys.argv[1:5]
    siti = [s.split(":") for s in sys.argv[5].split(",")]
    ril = Rilievo(ril_p)
    file = sorted(glob.glob(os.path.join(e36, "ecmwf_*.grb2")))
    righe = {s[0]: [] for s in siti}
    sel = {}
    for k, f in enumerate(file):
        run = dt.datetime.strptime(os.path.basename(f)[6:16], "%Y%m%d%H").replace(tzinfo=dt.timezone.utc)
        g = os.path.join(e912, os.path.basename(f).replace("ecmwf_", "e912_"))
        for ff in (f, g):
            if not os.path.exists(ff):
                continue
            try:
                lat, lon, C = leggi(ff)
            except Exception as ex:
                print("non letto", ff, ex); continue
            for nome, la0, lo0 in siti:
                if nome not in sel:
                    sel[nome] = nodi_largo(lat, lon, ril, float(la0), float(lo0))
                    print(nome, "nodi al largo:", int(sel[nome][0].sum()), "soglia km", sel[nome][1], flush=True)
                m = sel[nome][0]
                for (ch, st), U in C.items():
                    if ch != "10u" or ("10v", st) not in C:
                        continue
                    V = C[("10v", st)]
                    u, v = U[m], V[m]
                    S0 = float(np.median(np.hypot(u, v))); D0 = (np.degrees(np.arctan2(-u.mean(), -v.mean())) + 360) % 360
                    t = int((run + dt.timedelta(hours=int(st))).timestamp())
                    righe[nome].append((t, S0, D0, int(st)))
        if k % 100 == 0:
            print(k, len(file), flush=True)
    out = {}
    for nome, r in righe.items():
        r = np.array(sorted(set(r)))
        out[f"{nome}_t"] = r[:, 0].astype(np.int64); out[f"{nome}_S0"] = r[:, 1]; out[f"{nome}_D0"] = r[:, 2]; out[f"{nome}_passo"] = r[:, 3]
    np.savez_compressed(usc, **out)
    print("scritto", usc, {k: len(v) for k, v in righe.items()})


if __name__ == "__main__":
    main()
