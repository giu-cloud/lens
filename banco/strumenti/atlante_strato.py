"""Atlante dello strato marino per le aree del banco: per 8 direzioni del vento, la risposta dello strato sul rilievo fine
e su quello liscio del modello globale, portata sulla griglia del radar.

  python3 atlante_strato.py <rilievo .rgr> <cartella del metro (campi_*.npz)> <uscita> [H Fr passo_km liscio_km] [aree...]

Per ogni area scrive strato_<area>.npz con, per direzione (0, 45, ... 315 gradi di provenienza):
  loc     velocità dello strato nella cella ÷ vento imposto (rilievo fine)
  lett    velocità letta dai punti di mare a 0,25° sullo strato con il rilievo liscio ÷ vento imposto
  lettf   come lett, ma sullo strato con il rilievo fine
  spess   spessore dello strato nella cella ÷ H (rilievo fine)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rilievo import Rilievo
from strato_marino import al_punto, griglia, lettura_mare, risolvi

AREE = ["egadi", "eolie", "palermo", "messina", "malta", "napoli", "cagliari", "maddalena-bonifacio", "toscano", "genova"]


def main():
    ril_p, cart, usc = sys.argv[1:4]
    resto = sys.argv[4:]
    num = [x for x in resto if x.replace(".", "").isdigit()]
    aree = [x for x in resto if x in AREE] or AREE
    H, Fr, passo, liscio = (float(x) for x in (num + ["700", "1.0", "1.2", "4"][len(num):]))
    os.makedirs(usc, exist_ok=True)
    ril = Rilievo(ril_p)
    W = 8.0
    for area in aree:
        dst = os.path.join(usc, f"strato_{area}.npz")
        if os.path.exists(dst):
            continue
        t0 = time.time()
        z = np.load(os.path.join(cart, f"campi_{area}.npz"))
        lat, lon = z["lat"], z["lon"]
        g = griglia(ril, lat.min(), lat.max(), lon.min(), lon.max(), passo_km=passo, liscio_km=liscio)
        out = {k: np.zeros((8, len(lat), len(lon)), dtype="f4") for k in ("loc", "lett", "lettf", "spess")}
        for k in range(8):
            a = np.radians(45.0 * k)
            U, V = -W * np.sin(a), -W * np.cos(a)
            rf = risolvi(g, U, V, H, Fr)
            rl = risolvi(g, U, V, H, Fr, liscio=True)
            out["loc"][k] = al_punto(g, rf["S"], lat, lon) / W
            out["lett"][k] = lettura_mare(g, rl["S"], lat, lon) / W
            out["lettf"][k] = lettura_mare(g, rf["S"], lat, lon) / W
            out["spess"][k] = al_punto(g, rf["h"], lat, lon) / H
        np.savez(dst, H=H, Fr=Fr, passo=passo, liscio=liscio, **out)
        print(f"{area}: griglia {g['b'].shape}, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
