"""Atlante di Lens: la risposta del rilievo di un'area a venti di prova, calcolata una volta sola (metodo in METODO.md).

Con la stabilita' fissata, l'aggiustamento a conservazione della massa e' lineare nel campo di partenza: la risposta a
una somma di venti e' la somma delle risposte. Quindi per ogni classe di stabilita' (numero di Froude) si calcola come
il rilievo trasforma un vento da est di 1 m/s e uno da nord di 1 m/s (profilo logaritmico fino a 100 m, poi costante).
Sul tablet il vento fine di un'ora e' la combinazione delle due risposte, pesate con il vento del modello in quel punto,
tra le due classi di stabilita' piu' vicine; poi l'app aggiunge ridosso, fascia costiera e raffiche, che sono calcoli
leggeri. Cosi' Lens funziona in pochi secondi, senza rete, con qualsiasi GRIB (anche youmeteo).

Cosa si perde rispetto al calcolo completo: la stabilita' ora per ora dal profilo del modello (qui per classi) e il
vento in quota diverso da quello in basso (qui stessa direzione, profilo logaritmico).

Parametri: LENS_AREA "S,N,O,E", LENS_ATLANTE_ID (nome dell'area). Uscita: lens/uscita/atlante_<id>.lna (gzip).
Formato (little-endian): "LNA1"; int32 NX, NY, nx, ny, m, classi; float64 S, N, O, E del dominio e del riquadro;
float32 dx, H (quota di sbarramento), N efficace; float32 Fr di ogni classe; int16 quota (NY x NX, righe da nord);
float32 z0; uint8 mare; poi per ogni classe quattro campi int16 x 10000 sul riquadro (ny x nx, righe da nord):
u e v con vento da est, u e v con vento da nord (m/s per 1 m/s del modello a 10 m)."""
import gzip, os, struct, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LENS_PAGINA", "0")
from lens import calcola as C
from lens.solutore import livelli, aggiusta, al_suolo

# Classi di stabilita' (Fr dello strato dal mare alla cima dei rilievi). Sul tablet si interpola tra le due vicine.
CLASSI = [0.3, 0.5, 0.7, 1.0, 1.4]


def main():
    from scipy.ndimage import gaussian_filter
    A = C.A
    ident = "".join(c for c in (os.environ.get("LENS_ATLANTE_ID") or "area") if c.isalnum() or c == "-")[:40] or "area"
    t0 = time.time()
    print(f"atlante {ident}: dominio {A.S:.3f}-{A.N:.3f} N, {A.O:.3f}-{A.E:.3f} E · {A.NY} x {A.NX} celle da {A.dx:.0f} m")
    h, z0, mare, _, _ = C.terreno()
    terra = ~mare
    h_ril = float(max(200.0, min(1800.0, max(np.percentile(h[terra], 98), 0.85 * h.max()) if terra.any() else 200.0)))
    print(f"terreno: quota massima {h.max():.0f} m, sbarramento {h_ril:.0f} m, mare {mare.mean() * 100:.0f}%")
    dz, zc = livelli()
    m = A.m
    dentro = (slice(m, A.NY - m), slice(m, A.NX - m))
    ny, nx = A.NY - 2 * m, A.NX - 2 * m
    uno, zero = np.ones(h.shape), np.zeros(h.shape)
    risposte = []
    for fr in CLASSI:
        alfa2, hd, a_mid, a_alto = C.profilo_alfa2(zc, h_ril, fr, fr)
        memoria = {}
        campi = []
        for ue, ve in ((uno, zero), (zero, uno)):
            t1 = time.time()
            u3, v3 = C.iniziale(h, z0, zc, ue, ve, None, 0)
            # l'atlante si calcola una volta sola e si somma per classi sul tablet: qui la precisione piena (1e-5)
            u, v, w, aria, (d0, d1) = aggiusta(h, u3, v3, A.dx, dz, zc, alfa2=alfa2, tol=1e-5, memoria=memoria)
            del u3, v3
            us, vs = al_suolo(u, v, aria, h, zc, z0)
            del u, v, w, aria
            us, vs = gaussian_filter(us, 1.0), gaussian_filter(vs, 1.0)
            campi += [us[dentro], vs[dentro]]
            vel = np.hypot(us, vs)[dentro][mare[dentro]]
            print(f"Fr {fr}: bloccato fino a {hd:.0f} m, alfa2 {a_mid}/{a_alto} · vento da {'est' if ue is uno else 'nord'}: "
                  f"div {d0:.1e}->{d1:.1e} · mare {vel.min():.2f}/{np.median(vel):.2f}/{vel.max():.2f} · {time.time() - t1:.0f} s")
        risposte.append(campi)
    out = os.path.join(os.path.dirname(__file__), "uscita")
    os.makedirs(out, exist_ok=True)
    s_, n_, o_, e_ = A.interno
    buf = bytearray(b"LNA1")
    buf += struct.pack("<6i", A.NX, A.NY, nx, ny, m, len(CLASSI))
    buf += struct.pack("<8d", A.S, A.N, A.O, A.E, s_, n_, o_, e_)
    buf += struct.pack("<3f", A.dx, h_ril, C.N_MIN)
    buf += struct.pack(f"<{len(CLASSI)}f", *CLASSI)
    buf += np.clip(np.round(h), -32000, 32000).astype("<i2").tobytes()
    buf += z0.astype("<f4").tobytes()
    buf += mare.astype(np.uint8).tobytes()
    for campi in risposte:
        for c in campi:
            buf += np.clip(np.round(c * 10000.0), -32767, 32767).astype("<i2").tobytes()
    nome = os.path.join(out, f"atlante_{ident}.lna")
    with open(nome, "wb") as f:
        f.write(gzip.compress(bytes(buf), 6))
    print(f"atlante scritto: {os.path.getsize(nome) // 1024} KB, {len(CLASSI)} classi, {(time.time() - t0) / 60:.1f} minuti")


if __name__ == "__main__":
    main()
