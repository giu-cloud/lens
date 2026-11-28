"""La direzione del vento dalle strisce in un'immagine radar (metodo dei gradienti locali, Koch 2004, in forma semplice).

L'immagine è l'ampiezza a celle di 40 m (righe = azimut, colonne = distanza), con NaN su terra e fuori immagine.
  direzioni(img, passo_m=40, cella_m=200, finestra_m=2500) -> dict con, per finestra:
     r, c       il centro, in celle dell'immagine a 40 m
     ang        l'orientamento delle strisce nel piano dell'immagine, in gradi da 0 a 180,
                misurato dall'asse delle colonne verso l'asse delle righe
     q          la qualità, da 0 a 1 (quanto l'orientamento dei gradienti è netto)
     n          la quota di celle buone nella finestra
"""
import numpy as np
from scipy import ndimage as ndi


def liscia(a, sigma):
    """Media gaussiana che salta i NaN (convoluzione normalizzata)."""
    ok = np.isfinite(a)
    num = ndi.gaussian_filter(np.where(ok, a, 0.0), sigma, mode="nearest")
    den = ndi.gaussian_filter(ok.astype(float), sigma, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.3, num / den, np.nan)


def riduci(a, f):
    """Media a blocchi f x f che salta i NaN; NaN se meno di metà del blocco è buona."""
    ny, nx = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    b = a[:ny, :nx].reshape(ny // f, f, nx // f, f)
    ok = np.isfinite(b)
    n = ok.sum((1, 3)); s = np.where(ok, b, 0.0).sum((1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n >= f * f / 2, s / n, np.nan)


def prepara(img, passo_m=40.0, cella_m=200.0, fondo_m=4000.0, navi=4.0):
    """Toglie i punti troppo luminosi (navi, piattaforme), divide per il fondo lento (la luminosità che cambia con la
    distanza dal satellite e con la forza del vento) e porta l'immagine a celle di cella_m."""
    a = img.astype(float)
    f = max(1, int(round(cella_m / passo_m)))
    fondo = liscia(riduci(a, f), fondo_m / (passo_m * f))
    fondo_g = ndi.zoom(np.nan_to_num(fondo, nan=np.nanmedian(fondo)), f, order=1)[: a.shape[0], : a.shape[1]]
    if fondo_g.shape != a.shape:
        fondo_g = np.pad(fondo_g, ((0, a.shape[0] - fondo_g.shape[0]), (0, a.shape[1] - fondo_g.shape[1])), mode="edge")
    rel = a / fondo_g
    sd = np.nanstd(rel)
    forte = rel > 1 + navi * sd
    forte = ndi.binary_dilation(forte, iterations=3)
    rel[forte] = np.nan
    return riduci(rel, f), f


def direzioni(img, passo_m=40.0, cella_m=200.0, finestra_m=2500.0, passo_finestra_m=1000.0):
    b, f = prepara(img, passo_m, cella_m)
    b = liscia(b, 0.8)
    # gradienti (Scharr), solo dove tutto l'intorno 3x3 è buono
    gx = ndi.correlate(np.nan_to_num(b), np.array([[-3, 0, 3], [-10, 0, 10], [-3, 0, 3]]) / 32.0, mode="nearest")
    gy = ndi.correlate(np.nan_to_num(b), np.array([[-3, -10, -3], [0, 0, 0], [3, 10, 3]]) / 32.0, mode="nearest")
    ok = ndi.uniform_filter(np.isfinite(b).astype(float), 3, mode="constant") > 0.999
    g = (gx + 1j * gy) * ok
    m = np.abs(g)
    tetto = np.percentile(m[ok], 99) if ok.any() else 1.0   # i gradienti enormi sono bordi, non strisce
    ok &= m <= tetto
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(ok & (m > 0), g * g / m, 0.0)          # angolo doppio, peso = modulo
    w = max(3, int(round(finestra_m / (passo_m * f))))
    sz = ndi.uniform_filter(z.real, w, mode="constant") + 1j * ndi.uniform_filter(z.imag, w, mode="constant")
    sm = ndi.uniform_filter(np.abs(z), w, mode="constant")
    sn = ndi.uniform_filter(ok.astype(float), w, mode="constant")
    p = max(1, int(round(passo_finestra_m / (passo_m * f))))
    R, C = np.meshgrid(np.arange(p // 2, b.shape[0], p), np.arange(p // 2, b.shape[1], p), indexing="ij")
    Z = sz[R, C]; M = sm[R, C]; N = sn[R, C]
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(M > 0, np.abs(Z) / M, 0.0)
    grad = 0.5 * np.degrees(np.angle(Z))                    # orientamento dei gradienti, modulo 180
    ang = (grad + 90.0) % 180.0                             # le strisce sono perpendicolari ai gradienti
    return dict(r=(R + 0.5) * f, c=(C + 0.5) * f, ang=ang, q=q, n=N)


def finta(ny=1500, nx=1500, ang=30.0, passo_m=40.0, lunghezza_m=1000.0, forza=0.08, onde_m=220.0, onde_ang=110.0, onde=0.10,
          rumore=0.12, seme=0):
    """Un'immagine finta: strisce di direzione nota (larghezza irregolare), onde del mare più corte in un'altra direzione,
    luminosità che cala con la distanza, rumore moltiplicativo come quello del radar a 40 m."""
    rng = np.random.default_rng(seme)
    r, c = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    def lungo(a):   # coordinata perpendicolare a una direzione, in metri
        t = np.radians(a); return (-np.sin(t) * c + np.cos(t) * r) * passo_m
    # strisce: rumore filtrato lungo la direzione, non un seno perfetto
    base = rng.normal(size=(ny, nx))
    t = np.radians(ang)
    k = ndi.gaussian_filter(base, 2.0)
    k = ndi.rotate(ndi.gaussian_filter(ndi.rotate(k, ang, reshape=False, mode="reflect"), (lunghezza_m / passo_m / 6.0, 60.0), mode="reflect"), -ang, reshape=False, mode="reflect")
    k = k / k.std()
    img = 1.0 + forza * k + onde * np.sin(2 * np.pi * lungo(onde_ang) / onde_m)
    img *= np.linspace(1.3, 0.8, nx)[None, :]
    img *= 1.0 + rumore * rng.normal(size=(ny, nx))
    return np.clip(img, 0.05, None)


if __name__ == "__main__":
    # cancello 0: sulle immagini finte la direzione ritrovata deve stare entro 5 gradi
    # Il cancello vale per le celle di 200 m, che sono il metodo. Le celle di 100 m sono qui per mostrare perché no:
    # non tolgono le onde del mare (qui 220 m), che tirano la direzione dalla loro parte.
    peggio = 0.0
    for cella in (100.0, 200.0):
        for k, a in enumerate((0.0, 20.0, 45.0, 75.0, 90.0, 120.0, 160.0)):
            img = finta(ang=a, seme=k)
            img[:200, :300] = np.nan                         # un pezzo di terra
            d = direzioni(img, cella_m=cella)
            buone = (d["n"] > 0.8)
            e = (d["ang"][buone] - a + 90) % 180 - 90
            med = np.degrees(0.5 * np.angle(np.mean(np.exp(2j * np.radians(e)))))
            print(f"celle {cella:.0f} m, strisce a {a:5.1f}: errore medio {med:+5.1f}, mediano assoluto {np.median(np.abs(e)):4.1f}, entro 10 gradi {np.mean(np.abs(e) < 10) * 100:3.0f}%, qualità mediana {np.median(d['q'][buone]):.2f}")
            if cella == 200.0:
                peggio = max(peggio, abs(med), np.median(np.abs(e)))
    # senza strisce: la qualità deve essere bassa
    img = finta(forza=0.0, seme=9); d = direzioni(img)
    print(f"senza strisce: qualità mediana {np.median(d['q'][d['n'] > 0.8]):.2f}")
    print("CANCELLO 0:", "passa" if peggio <= 5 else "NON passa", f"(peggiore {peggio:.1f})")
