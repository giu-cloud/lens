"""Prova indipendente dalla fase della marea del modello: residuo del radar su marea all'ora del passaggio e 3 ore prima (un quarto di ciclo)."""

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
k3 = k - 3
cq = tu[k3, jj, ii] * wx + tv[k3, jj, ii] * wy
T = np.load('/home/claude/lavoro/tabella6.npz', mmap_mode='r'); sel = np.asarray(T['area']) == 3
vero = np.load('/home/claude/lavoro/radar_vero_tabella6.npy', mmap_mode='r')[sel]
ok = okt & (dkm < 3) & np.isfinite(cwt) & np.isfinite(cq) & np.isfinite(vero) & (vero > 1) & (mod > 3)
y = np.log(vero / mod); key = (j * 100 + i) * 3 + reg
_, inv = np.unique(key[ok], return_inverse=True); m_ = np.bincount(inv, y[ok]) / np.bincount(inv)
res = np.full(len(y), np.nan); res[np.where(ok)[0]] = y[ok] - m_[inv]
for nome, m in (("tutto", ok), ("stretta", ok & (la > 38.15) & (lo > 15.55)), ("vento oltre 6 m/s", ok & (mod > 6))):
    X = np.c_[-cwt[m] / mod[m], -cq[m] / mod[m], -cwl[m] / mod[m], np.ones(m.sum())]; Y = res[m]
    b = np.linalg.lstsq(X, Y, rcond=None)[0]
    ps = sc[m]; up = np.unique(ps); idx = {p: np.where(ps == p)[0] for p in up}; rng = np.random.default_rng(2); bs = []
    for _ in range(300):
        q = np.concatenate([idx[p] for p in rng.choice(up, len(up))]); bs.append(np.linalg.lstsq(X[q], Y[q], rcond=None)[0][:2])
    bs = np.array(bs); amp = np.hypot(b[0], b[1]); ampb = np.hypot(bs[:, 0], bs[:, 1])
    print(f"  {nome:18s} in fase {b[0]:+.2f} ± {bs[:,0].std():.2f}, 3 h prima {b[1]:+.2f} ± {bs[:,1].std():.2f}; "
          f"ampiezza {amp:.2f} (90% dei ricampionamenti tra {np.percentile(ampb,5):.2f} e {np.percentile(ampb,95):.2f}); fase {np.degrees(np.arctan2(b[1], b[0])):+.0f}°")
