"""Lettore del rilievo di Shift (formato RGR2 di core/Relief.kt): quote in metri su più livelli, a riquadri compressi.

  r = Rilievo("europa-....rgr"); r.livelli; h = r.griglia(livello, sud, nord, ovest, est)  -> (lat, lon, quote)
Le righe vanno da sud a nord; 0 = mare.
"""
import struct
import zlib

import numpy as np


class Rilievo:
    def __init__(self, percorso):
        self.b = open(percorso, "rb").read()
        assert self.b[:4] == b"RGR2", "non è un RGR2"
        n = struct.unpack(">i", self.b[4:8])[0]
        pos = 8
        self.livelli = []
        for _ in range(n):
            lat0, lon0, cell = struct.unpack(">ddd", self.b[pos:pos + 24]); pos += 24
            rows, cols, tile, ty, tx = struct.unpack(">iiiii", self.b[pos:pos + 20]); pos += 20
            ol = np.frombuffer(self.b[pos:pos + 8 * ty * tx], dtype=">i4").reshape(-1, 2); pos += 8 * ty * tx
            self.livelli.append(dict(lat0=lat0, lon0=lon0, cell=cell, rows=rows, cols=cols, tile=tile, ty=ty, tx=tx,
                                     off=ol[:, 0].copy(), len=ol[:, 1].copy()))
        self.dati = pos

    def riquadro(self, lv, ty, tx):
        L = self.livelli[lv]; t = L["tile"]
        k = ty * L["tx"] + tx
        if L["len"][k] == 0:
            return np.zeros((t, t), dtype=np.int32)
        s = self.dati + int(L["off"][k])
        raw = np.frombuffer(zlib.decompress(self.b[s:s + int(L["len"][k])]), dtype=np.uint8)
        c = t * t
        d = raw[:c].astype(np.int64) | (raw[c:2 * c].astype(np.int64) << 8)
        return (np.cumsum(d.reshape(t, t), axis=1) & 0xFFFF).astype(np.uint16).astype(np.int16).astype(np.int32)

    def griglia(self, lv, s, n, w, e):
        L = self.livelli[lv]; c = L["cell"]; t = L["tile"]
        j0 = max(int(np.floor((s - L["lat0"]) / c)), 0); j1 = min(int(np.ceil((n - L["lat0"]) / c)), L["rows"])
        i0 = max(int(np.floor((w - L["lon0"]) / c)), 0); i1 = min(int(np.ceil((e - L["lon0"]) / c)), L["cols"])
        out = np.zeros((j1 - j0, i1 - i0), dtype=np.int32)
        for ty in range(j0 // t, (j1 - 1) // t + 1):
            for tx in range(i0 // t, (i1 - 1) // t + 1):
                if ty >= L["ty"] or tx >= L["tx"]:
                    continue
                q = self.riquadro(lv, ty, tx)
                ya, yb = max(j0, ty * t), min(j1, (ty + 1) * t); xa, xb = max(i0, tx * t), min(i1, (tx + 1) * t)
                out[ya - j0:yb - j0, xa - i0:xb - i0] = q[ya - ty * t:yb - ty * t, xa - tx * t:xb - tx * t]
        lat = L["lat0"] + (np.arange(j0, j1) + 0.5) * c
        lon = L["lon0"] + (np.arange(i0, i1) + 0.5) * c
        return lat, lon, out
