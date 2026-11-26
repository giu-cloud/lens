"""Le aree fini allargate (0.57.1) a confronto con Lens 0.56 sull'area grande, per un vento uniforme (15 nodi, Fr 0,7).
Le aree fini sono messe sopra l'area grande come le mette l'app (LensWind una sopra l'altra: raccordo di 1,5 miglia
sui bordi di ognuna). Campi da banco/cfd/Campi.kt: <pref_fini><id>.bin e <pref_madri><area>.bin, due venti (315°, 270°).
Uso: python3 figura_aree_fini.py <cartella campi> <pref_fini> <pref_madri> <atlanti fini> <atlanti grandi> <png> <aree,...> <venti: 0,1>"""
import json, numpy as np, sys, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
qui = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [qui, os.path.join(qui, '..', 'cfd'), os.path.join(qui, '..', '..', 'lens')]
from leggi_campi import leggi
from controlla_atlante import leggi as leggi_lna
import aree
C, PF, PM, AF, AG, png, madri, venti = sys.argv[1:9]
madri = madri.split(','); venti = [int(x) for x in venti.split(',')]
NOMI_VENTI = ["da 315°", "da 270°"]
# VENTI=<json {area: [da1, da2]}>: le direzioni con cui sono stati calcolati i campi di quell'area (default 315 e 270)
VENTI = json.load(open(os.environ["VENTI"])) if os.environ.get("VENTI") else {}
RACCORDO_KM = 1.5 * 1.852
def mare_di(p):
    A = leggi_lna(p); m = A['m']; return A['mare'][m:A['NY'] - m, m:A['NX'] - m]
def interp(la, lo, Z, LA, LO, metodo='linear'):
    f = RegularGridInterpolator((la[::-1], lo), Z[::-1].astype(float), method=metodo, bounds_error=False, fill_value=np.nan)
    return f(np.stack([LA.ravel(), LO.ravel()], -1)).reshape(LA.shape)
def mosaico(p, r):
    """Il campo che l'app mostra: l'area grande (primo Lens) con sopra le aree fini, raccordate sui bordi."""
    S, N, O, E = aree.AREE[p]; kml = 111.19 * np.cos(np.radians((S + N) / 2))
    la = np.linspace(N, S, int((N - S) * 111.19 / 0.1)); lo = np.linspace(O, E, int((E - O) * kml / 0.1))
    LA, LO = np.meshgrid(la, lo, indexing='ij')
    gl, go, F = leggi(f"{C}/{PM}{p}.bin", 2); (u, v) = F[r][1]
    U, V = interp(gl, go, u, LA, LO), interp(gl, go, v, LA, LO)
    M = interp(gl, go, mare_di(f"{AG}/atlante_{p}.lna"), LA, LO, 'nearest') > 0.5
    for k in [k for k, q in aree.PADRE.items() if q == p]:
        fl, fo, Ff = leggi(f"{C}/{PF}{k}.bin", 2); (uf, vf) = Ff[r][1]
        s, n, o, e = aree.FINI[k]
        bordo = np.minimum(np.minimum(LA - s, n - LA) * 111.19, np.minimum(LO - o, e - LO) * kml)
        w = np.clip(bordo / RACCORDO_KM, 0, 1)
        dentro = bordo > 0
        Uf, Vf = interp(fl, fo, uf, LA, LO), interp(fl, fo, vf, LA, LO)
        Mf = interp(fl, fo, mare_di(f"{AF}/atlante_{k}.lna"), LA, LO, 'nearest') > 0.5
        ok = dentro & np.isfinite(Uf)
        # come LensWind: la forza si raccorda a parte, la direzione sulle componenti
        sf, sb = np.hypot(Uf, Vf), np.hypot(U, V)
        uu, vv = w * Uf + (1 - w) * U, w * Vf + (1 - w) * V
        ss = w * sf + (1 - w) * sb; k_ = ss / (np.hypot(uu, vv) + 1e-9)
        U = np.where(ok, uu * k_, U); V = np.where(ok, vv * k_, V); M = np.where(ok, Mf, M)
    return la, lo, U, V, M
def pannello(a, la, lo, U, V, M, passo_km):
    S = np.hypot(U, V) * 1.944
    im = a.pcolormesh(lo, la, np.where(M, S, np.nan), cmap='turbo', vmin=0, vmax=24, shading='auto', rasterized=True)
    a.contour(lo, la, M.astype(float), [0.5], colors='k', linewidths=0.5)
    k = max(1, int(round(passo_km / (abs(la[0] - la[1]) * 111.19))))
    LA, LO = np.meshgrid(la, lo, indexing='ij'); sel = np.zeros_like(M); sel[::k, ::k] = True; sel &= M & np.isfinite(U)
    sp = np.hypot(U, V) + 1e-9
    a.quiver(LO[sel], LA[sel], (U / sp)[sel], (V / sp)[sel], scale=70, width=0.0012, color='k', alpha=0.7)
    return im
def riquadro(a, b, **kw):
    S, N, O, E = b; a.plot([O, E, E, O, O], [S, S, N, N, S], **kw)
righe = [(p, r) for p in madri for r in venti]
fig, ax = plt.subplots(len(righe), 2, figsize=(22, 8.5 * len(righe)), squeeze=False)
for i, (p, r) in enumerate(righe):
    gl, go, F = leggi(f"{C}/{PM}{p}.bin", 2); (u0, v0) = F[r][0] if os.environ.get("SINISTRA") == "0.56" else F[r][1]
    S_, N_, O_, E_ = aree.AREE[p]; passo = max(1.2, (E_ - O_) * 111.19 * np.cos(np.radians((S_ + N_) / 2)) / 40)
    im = pannello(ax[i, 0], gl, go, u0, v0, mare_di(f"{AG}/atlante_{p}.lna"), passo)
    la, lo, U, V, M = mosaico(p, r)
    pannello(ax[i, 1], la, lo, U, V, M, passo)
    for k in [k for k, q in aree.PADRE.items() if q == p]:
        riquadro(ax[i, 1], aree.FINI[k], color='w', lw=0.8, ls=':')
    S, N, O, E = aree.AREE[p]
    ax[i, 0].set_title(f"{p}, {('da %g°' % VENTI[p][r]) if p in VENTI else NOMI_VENTI[r]}, 15 nodi inventati, aria stabile — Lens sull'area grande, maglia grossa" if os.environ.get("SINISTRA") != "0.56" else "15 nodi, Fr 0,7 — Lens 0.56 sull'area grande", fontsize=10)
    ax[i, 1].set_title(f"{p} — lo stesso Lens con le aree fini (puntini bianchi), raccordate come nell'app", fontsize=10)
    for a in ax[i]:
        a.set_xlim(O, E); a.set_ylim(S, N); a.set_aspect(1 / np.cos(np.radians((S + N) / 2)))
plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.4, label='nodi', pad=0.01)
plt.savefig(png, dpi=55, bbox_inches='tight')
