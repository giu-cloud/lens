"""Lo Stretto con la corrente (26/11): il vento sul fondo (C ora, livello su 5 km), la corrente del modello Copernicus (con
le maree) alla stessa ora, il vento sull'acqua (vento meno corrente: quello che sente la barca, come Router.waterWind di
Shift) e la differenza. Ora vera, previsione ECMWF. Uso: python3 figura_corrente.py <corrente.bin> <correnti.npz> <png>"""
import sys, datetime as dt, numpy as np, matplotlib
from zoneinfo import ZoneInfo
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator
sys.argv, argv = sys.argv[:1] + ['-', '0', '-'], sys.argv
exec(open(__file__.replace('figura_corrente.py', 'figura_schermate.py')).read().replace("\nmain()\n", "\n"))
CORR = LinearSegmentedColormap.from_list('corr', ['#E8F1F8', '#6BAED6', '#08519C', '#08306B'])
box, A = leggi(argv[1]); f = A[0]; S, N, O, E = box; ny, nx = f.shape[1:]
lat = N - (np.arange(ny) + 0.5) * (N - S) / ny; lon = O + (np.arange(nx) + 0.5) * (E - O) / nx
LA, LO = np.meshgrid(lat, lon, indexing='ij'); mare = f[8] > 0.5
Z = np.load(argv[2]); tc = Z['t']
t = int(open('/tmp/sch/ora_corrente.txt').read())
k = int(np.argmin(np.abs(tc - t))); assert abs(tc[k] - t) <= 1800
def riempi(a):
    a = a.copy(); ok = np.isfinite(a); jj, ii = np.nonzero(ok)
    for j, i in zip(*np.nonzero(~ok)):
        q = np.argmin((jj - j) ** 2 + (ii - i) ** 2); a[j, i] = a[jj[q], ii[q]]
    return a
cu, cv = riempi(Z['u'][k].astype(float)), riempi(Z['v'][k].astype(float))
fi = lambda a: RegularGridInterpolator((Z['lat'], Z['lon']), a, bounds_error=False, fill_value=None)(np.stack([LA.ravel(), LO.ravel()], -1)).reshape(LA.shape)
Cu, Cv = fi(cu), fi(cv)
sC = np.hypot(f[3], f[4]); k5 = f[7]
Wu, Wv = f[3] * k5, f[4] * k5                  # C ora: vento sul fondo
Au, Av = Wu - Cu, Wv - Cv                      # vento sull'acqua
pann = [("Vento sul fondo: C ora (Lens)", np.hypot(Wu, Wv), Wu, Wv, VENTO, 0, 32 / KN),
        ("Corrente del modello (Copernicus, con le maree)", np.hypot(Cu, Cv), Cu, Cv, CORR, 0, 4 / KN),
        ("Vento sull'acqua: quello che sente la barca\n(vento meno corrente)", np.hypot(Au, Av), Au, Av, VENTO, 0, 32 / KN),
        ("Differenza: sull'acqua meno sul fondo", np.hypot(Au, Av) - np.hypot(Wu, Wv), None, None, DIFF, -4 / KN, 4 / KN)]
d = dt.datetime.fromtimestamp(t, dt.timezone.utc); loc = d.astimezone(ZoneInfo('Europe/Rome'))
fig, ax = plt.subplots(1, 4, figsize=(21, 9.6))
for a, (nome, s, u, v, cm, lo_, hi) in zip(ax, pann):
    a.set_facecolor(TERRA); a.set_xlim(O, E); a.set_ylim(S, N); a.set_aspect(1 / np.cos(np.radians((S + N) / 2))); a.set_xticks([]); a.set_yticks([])
    im = a.pcolormesh(lon, lat, np.where(mare, s * KN, np.nan), cmap=cm, vmin=lo_ * KN, vmax=hi * KN, shading='auto', rasterized=True)
    a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
    if u is not None:
        p = 15; sel = np.zeros_like(mare); sel[p // 2::p, p // 2::p] = True; sel &= mare & np.isfinite(u)
        sp = np.hypot(u, v) + 1e-9
        a.quiver(LO[sel], LA[sel], (u / sp)[sel], (v / sp)[sel], scale=38, width=0.004, color='#1F2D3A', alpha=0.65)
    q = mare & (LA > 38.215) & (LA < 38.265) & (LO > 15.58) & (LO < 15.66)
    a.set_title(nome + (f"\nnella stretta, in media: {np.nanmean(s[q]) * KN:+.1f} nodi; al massimo {np.nanmax(s[q]) * KN:+.1f}" if lo_ < 0 else f"\nnella stretta, in media: {np.nanmean(s[q]) * KN:.1f} nodi").replace('.', ','), fontsize=12)
    cb = fig.colorbar(im, ax=a, shrink=0.6, pad=0.02, orientation='horizontal'); cb.set_label('nodi')
    for la_, lo2, nm in ((38.19, 15.555, "Messina"), (38.11, 15.65, "Reggio"), (38.26, 15.61, ""), (38.22, 15.637, "")):
        a.plot(lo2, la_, 'o', ms=2.5, color='#1F2D3A')
        if nm: a.text(lo2 + (0.012 if nm == 'Reggio' else -0.008), la_, nm, fontsize=8, ha='left' if nm == 'Reggio' else 'right', color='#1F2D3A')
fig.suptitle(f"Stretto di Messina, {loc:%d/%m/%Y %H:%M} ora italiana: maestrale di circa 20 nodi e montante di 3,4 nodi alla stretta (corrente contro il vento).\n"
             f"Previsione ECMWF del 17/02 12 UTC; corrente Copernicus della stessa ora. Esempio scelto tra le ore con corrente forte e vento forte.", fontsize=13, y=1.0)
fig.savefig(argv[3], dpi=70, bbox_inches='tight')
print(argv[3])
