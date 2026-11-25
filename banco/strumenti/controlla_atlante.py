"""Controllo di un atlante di Lens (LNA1): intestazione, terreno, e la risposta sul mare ai venti di prova.

  python3 controlla_atlante.py atlante_<area>.lna [altro.lna ...]

Per ogni classe di stabilità e per i due venti di prova (u = 1 m/s e v = 1 m/s, come in lens/atlante.py): la forza sul mare (mediana, 5° e 95° percentile, in
rapporto al vento del modello) e la rotazione media. Valori plausibili: mediana vicina a 1, poche celle sotto 0,5 o
sopra 1,5, rotazioni di pochi gradi al largo.
"""
import gzip
import struct
import sys

import numpy as np


def leggi(p):
    b = gzip.open(p).read(); o = 4
    assert b[:4] == b"LNA1", "non è LNA1"
    NX, NY, nx, ny, m, k = struct.unpack_from("<6i", b, o); o += 24
    S, N, O, E, s, n, w, e = struct.unpack_from("<8d", b, o); o += 64
    dx, H, Neff = struct.unpack_from("<3f", b, o); o += 12
    fr = struct.unpack_from(f"<{k}f", b, o); o += 4 * k
    h = np.frombuffer(b, "<i2", NX * NY, o).reshape(NY, NX); o += 2 * NX * NY
    z0 = np.frombuffer(b, "<f4", NX * NY, o).reshape(NY, NX); o += 4 * NX * NY
    mare = np.frombuffer(b, "u1", NX * NY, o).reshape(NY, NX).astype(bool); o += NX * NY
    campi = []
    for _ in range(k):
        c = []
        for _ in range(4):
            c.append(np.frombuffer(b, "<i2", nx * ny, o).reshape(ny, nx) / 10000.0); o += 2 * nx * ny
        campi.append(c)
    assert o == len(b), f"lunghezza: letti {o} byte su {len(b)}"
    return dict(NX=NX, NY=NY, nx=nx, ny=ny, m=m, box=(s, n, w, e), dominio=(S, N, O, E), dx=dx, H=H, fr=fr, h=h, mare=mare, campi=campi)


def main():
    for p in sys.argv[1:]:
        A = leggi(p)
        m = A["m"]; mare = A["mare"][m:A["NY"] - m, m:A["NX"] - m]
        print(f"{p}: riquadro {A['box']}, dominio {A['NX']}x{A['NY']} celle da {A['dx']:.0f} m, margine {m}, "
              f"quota massima {A['h'].max()} m, sbarramento {A['H']:.0f} m, mare {mare.mean() * 100:.0f}% del riquadro")
        for fr, (ue, ve, un, vn) in zip(A["fr"], A["campi"]):
            righe = []
            for nome, u, v, (u0, v0) in (("u=1", ue, ve, (1, 0)), ("v=1", un, vn, (0, 1))):
                vel = np.hypot(u, v)[mare]
                rot = np.degrees(np.arctan2(u * v0 - v * u0, u * u0 + v * v0))[mare]
                righe.append(f"{nome}: forza {np.median(vel):.2f} ({np.percentile(vel, 5):.2f}-{np.percentile(vel, 95):.2f}), "
                             f"rotazione media {np.mean(rot):+.1f}°")
            print(f"  Fr {fr:.1f}  " + "  |  ".join(righe))


if __name__ == "__main__":
    main()
