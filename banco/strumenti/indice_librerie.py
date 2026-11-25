"""Rifà indice.json delle librerie LRA1 dalla cartella: file, area, modello, celle, passaggi, byte, sha256, riquadro, attiva.

  python3 indice_librerie.py <cartella librerie> [--costruite AAAA-MM-GG] [--radar-fino-a AAAA-MM-GG]

Tiene i campi generali dell'indice che c'è già (pesi dell'insieme, nota) e il campo "attiva" delle librerie già presenti
(false = spenta da remoto: l'app non la applica).
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from libreria_radar import Libreria


def main():
    cart = sys.argv[1]
    arg = dict(zip(sys.argv[2::2], sys.argv[3::2]))
    p = os.path.join(cart, "indice.json")
    vecchio = json.load(open(p)) if os.path.exists(p) else {}
    attive = {x["file"]: x.get("attiva", True) for x in vecchio.get("librerie", [])}
    lib = []
    for f in sorted(os.listdir(cart)):
        if not f.endswith(".lra"):
            continue
        q = os.path.join(cart, f)
        L = Libreria(q)
        la = L.lat0 + (L.celle // L.nx) * L.dlat; lo = L.lon0 + (L.celle % L.nx) * L.dlon
        lib.append(dict(file=f, area=L.area, modello=L.modello, celle=int(L.nc), passaggi=int(L.ns), byte=os.path.getsize(q),
                        sha256=hashlib.sha256(open(q, "rb").read()).hexdigest(),
                        riquadro=[round(float(la.min()), 4), round(float(la.max()), 4), round(float(lo.min()), 4), round(float(lo.max()), 4)],
                        attiva=attive.get(f, True)))
    out = {k: v for k, v in vecchio.items() if k != "librerie"}
    out.setdefault("formato", "LRA1")
    if "--costruite" in arg:
        out["costruite"] = arg["--costruite"]
    if "--radar-fino-a" in arg:
        out["radar_fino_a"] = arg["--radar-fino-a"]
    out["librerie"] = lib
    json.dump(out, open(p, "w"), indent=1, ensure_ascii=False)
    print(f"{len(lib)} librerie, {sum(x['byte'] for x in lib) / 1e6:.2f} MB -> {p}")


if __name__ == "__main__":
    main()
