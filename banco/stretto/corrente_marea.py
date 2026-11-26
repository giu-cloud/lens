"""Come corrente_radar.py, ma la corrente è divisa in marea (astronomica, non dipende dal vento) e resto (filtro di Godin 24-24-25 h).
Solo la marea dà una prova pulita: il resto è spinto dal vento stesso e confonde."""
import numpy as np
Z = np.load('/home/claude/dati/correnti/correnti_stretto.npz')
tc, cla, clo, cu, cv = Z['t'], Z['lat'], Z['lon'], Z['u'].astype(float), Z['v'].astype(float)
assert np.all(np.diff(tc) == 3600), "ore non consecutive"
def media(a, n):
    k = np.ones(n) / n
    f = lambda s: np.convolve(np.nan_to_num(s), k, mode='same')
    return np.apply_along_axis(f, 0, a)
def godin(a): return media(media(media(a, 24), 24), 25)
lu, lv = godin(cu), godin(cv); tu, tv = cu - lu, cv - lv
mare = np.isfinite(cu[0])
s = np.hypot(tu, tv)
print("marea alla stretta: velocità 95%:", round(float(np.nanpercentile(s[:, 7, 5], 95)) * 1.944, 2), "nodi; resto 95%:",
      round(float(np.nanpercentile(np.hypot(lu, lv)[:, 7, 5], 95)) * 1.944, 2), "nodi")
T = np.load('/home/claude/lavoro/tabella6.npz', mmap_mode='r'); sel = np.asarray(T['area']) == 3
g = lambda k: np.asarray(T[k][sel])
vero = np.load('/home/claude/lavoro/radar_vero_tabella6.npy', mmap_mode='r')[sel]
C = np.load('/home/claude/lavoro/v_scia/campi_messina.npz'); rlat, rlon = C['lat'], C['lon']
j, i = g('j').astype(int), g('i').astype(int); mod, dirc, D0, t, dist, sc = g('mod'), g('dir'), g('D0'), g('t'), g('dist'), g('scena')
la, lo = rlat[j], rlon[i]
k = np.clip(np.searchsorted(tc, t), 1, len(tc) - 1); k = np.where(np.abs(tc[k - 1] - t) < np.abs(tc[k] - t), k - 1, k)
okt = (np.abs(tc[k] - t) <= 1800) & (k > 40) & (k < len(tc) - 40)
jj = np.abs(cla[None, :] - la[:, None]).argmin(1); ii = np.abs(clo[None, :] - lo[:, None]).argmin(1)
dkm = np.hypot((cla[jj] - la) * 111.19, (clo[ii] - lo) * 111.19 * np.cos(np.radians(38.1)))
a = np.radians(dirc); wx, wy = -np.sin(a), -np.cos(a)
cwt = tu[k, jj, ii] * wx + tv[k, jj, ii] * wy; cwl = lu[k, jj, ii] * wx + lv[k, jj, ii] * wy
reg = np.where((D0 >= 315) | (D0 < 45), 0, np.where((D0 >= 135) & (D0 < 225), 1, 2))
ok = okt & (dkm < 3) & np.isfinite(cwt) & np.isfinite(vero) & (vero > 1) & (mod > 3)
y = np.log(vero / mod)
key = (j * 100 + i) * 3 + reg
_, inv = np.unique(key[ok], return_inverse=True)
m_ = np.bincount(inv, y[ok]) / np.bincount(inv)
res = np.full(len(y), np.nan); res[np.where(ok)[0]] = y[ok] - m_[inv]
xt, xl = -cwt / mod, -cwl / mod
def blocchi(m, nome):
    # due variabili insieme; errore della pendenza con ricampionamento per passaggio (le righe di un passaggio non sono indipendenti)
    X = np.c_[xt[m], xl[m], np.ones(m.sum())]; Y = res[m]
    b = np.linalg.lstsq(X, Y, rcond=None)[0]
    ps = sc[m]; up = np.unique(ps); rng = np.random.default_rng(1); bs = []
    idx = {p: np.where(ps == p)[0] for p in up}
    for _ in range(300):
        q = np.concatenate([idx[p] for p in rng.choice(up, len(up))])
        bs.append(np.linalg.lstsq(X[q], Y[q], rcond=None)[0][:2])
    e = np.std(bs, 0)
    ct = cwt[m]; f = np.abs(ct) > 0.5144
    contro = res[m][f & (ct < 0)]; favore = res[m][f & (ct > 0)]
    pc = lambda r: f"{100*(np.exp(r.mean())-1):+.1f}% ({len(np.unique(ps[f][(ct[f]<0) if r is contro else (ct[f]>0)]))} pass.)" if len(r) > 50 else "  -"
    print(f"  {nome:24s} righe {m.sum():6d}  marea: pendenza {b[0]:+.2f} ± {e[0]:.2f}   resto: {b[1]:+.2f} ± {e[1]:.2f}   "
          f"marea oltre 1 nodo contro il vento {pc(contro)}, a favore {pc(favore)}")
print("pendenza 1 = il radar vede il vento sull'acqua (vento meno corrente); 0 = non la vede")
for nome, m in (("tutto", ok), ("vento da nord", ok & (reg == 0)), ("vento da sud", ok & (reg == 1)), ("entro 3 km", ok & (dist < 3)),
                ("stretta (oltre 38,15)", ok & (la > 38.15) & (lo > 15.55)), ("vento oltre 6 m/s", ok & (mod > 6)), ("vento 3-6 m/s", ok & (mod <= 6))):
    if m.sum() > 200: blocchi(m, nome)
np.savez_compressed('/home/claude/lavoro/stretto/marea_righe.npz', res=res[ok], xt=xt[ok], xl=xl[ok], cwt=cwt[ok], mod=mod[ok], reg=reg[ok], la=la[ok], lo=lo[ok], sc=sc[ok])
