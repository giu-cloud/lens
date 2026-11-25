"""Legge la pagina della carta di Pietro (pietro.html, release "pietro" di giu-cloud/pogoda): ore, forza, raffica e direzione.

  from leggi_pietro import leggi; P = leggi("pietro.html")
  P["ore"], P["nx"], P["ny"], P["mare"] (ny, nx), P["V"] (ore, ny, nx) nodi, P["G"] raffica, P["D"] direzione (da dove, gradi)

Come nella pagina: D.bin è gzip in base64, tre campi (forza e raffica in mezzi nodi, direzione in 128 passi), ognuno ora
dopo ora come differenze modulo 256; P (D.gz) ha le ore e la maschera del mare (un bit per cella).
"""
import base64
import gzip
import json
import re

import numpy as np


def leggi(p):
    t = open(p, encoding="utf-8").read()
    i = t.index("const D=") + len("const D=")
    D, _ = json.JSONDecoder().raw_decode(t[i:])
    P = json.loads(gzip.decompress(base64.b64decode(D["gz"])))
    B = np.frombuffer(gzip.decompress(base64.b64decode(D["bin"])), np.uint8)
    nx, ny = D["nx"], D["ny"]; N = nx * ny; NH = len(P["ore"])
    B = B[:3 * NH * N].reshape(3, NH, N)
    B = (np.cumsum(B.astype(np.int64), axis=1) % 256).astype(np.uint8)
    m = np.frombuffer(base64.b64decode(P["mare"]), np.uint8)
    mare = np.unpackbits(m)[:N].reshape(ny, nx).astype(bool)
    ore = [o["t"] for o in P["ore"]]
    return dict(D=D, ore=ore, nx=nx, ny=ny, mare=mare, interno=D.get("interno"),
                V=B[0].reshape(NH, ny, nx) / 2.0, G=B[1].reshape(NH, ny, nx) / 2.0, Dir=B[2].reshape(NH, ny, nx) * 360.0 / 128)
