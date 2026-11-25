"""Unisce gli elenchi delle stazioni scaricate (GHCNh, porti ISPRA, METAR di IEM, Météo-France) in uno solo.

  python3 stazioni_elenco.py <cartella stazioni> <uscita.json>

Le stesse stazioni compaiono in più fonti (per esempio Cap Pertusato in GHCNh e in Météo-France, Palermo Punta Raisi in
GHCNh e in IEM). Se due stazioni distano meno di 2,5 km si tiene quella della fonte migliore, nell'ordine:
Météo-France (ogni ora, media di 10 minuti a 10 m, controllata), porti ISPRA, GHCNh, IEM; a parità, quella con più misure.
Si tengono solo le stazioni con almeno 2000 misure.
"""
import json
import os
import sys

import numpy as np

ORDINE = {"mf": 0, "rmn": 1, "gh": 2, "iem": 3}


def main():
    cart, usc = sys.argv[1:3]
    tutte = []
    for fonte in ORDINE:
        p = os.path.join(cart, f"stazioni_{fonte}.json")
        if not os.path.exists(p):
            continue
        for s in json.load(open(p)):
            if s.get("misure", 0) >= 2000 and os.path.exists(os.path.join(cart, f"{fonte}_{s['codice']}.npz")):
                tutte.append(dict(s, fonte=fonte))
    tutte.sort(key=lambda s: (ORDINE[s["fonte"]], -s["misure"]))
    tenute = []
    for s in tutte:
        doppia = None
        for t in tenute:
            d = 111.2 * np.hypot(s["lat"] - t["lat"], (s["lon"] - t["lon"]) * np.cos(np.radians(t["lat"])))
            if d < 2.5:
                doppia = t; break
        if doppia:
            doppia.setdefault("doppioni", []).append(f"{s['fonte']}_{s['codice']}")
            continue
        tenute.append(s)
    json.dump(tenute, open(usc, "w"), indent=1, ensure_ascii=False)
    for a in sorted(set(s["area"] for s in tenute)):
        print(a, ", ".join(f"{s['fonte']}:{s['nome'][:22]}" for s in tenute if s["area"] == a))
    print(f"{len(tenute)} stazioni (da {len(tutte)})")


if __name__ == "__main__":
    main()
