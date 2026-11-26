"""Il radar nello Stretto vede il vento rispetto all'acqua? Residuo del radar contro la corrente lungo il vento (26/11).
Residuo = log(vero/modello) meno la sua media nella stessa cella e nello stesso regime (toglie il getto medio).
x = −(corrente lungo la direzione del vento) / vento del modello: se il radar misura il vento sull'acqua, residuo ≈ x."""
import numpy as np, sys
Z = np.load('/home/claude/dati/correnti/correnti_stretto.npz')
tc, cla, clo, cu, cv = Z['t'], Z['lat'], Z['lon'], Z['u'], Z['v']
T = np.load('/home/claude/lavoro/tabella6.npz', mmap_mode='r')
sel = np.asarray(T['area']) == 3
g = lambda k: np.asarray(T[k][sel])
vero = np.load('/home/claude/lavoro/radar_vero_tabella6.npy', mmap_mode='r')[sel]
C = np.load('/home/claude/lavoro/v_scia/campi_messina.npz'); rlat, rlon = C['lat'], C['lon']
j, i = g('j').astype(int), g('i').astype(int); mod, dirc, D0, t, dist, sc = g('mod'), g('dir'), g('D0'), g('t'), g('dist'), g('scena')
la, lo = rlat[j], rlon[i]
# corrente all'ora più vicina, punto del modello più vicino (entro 3 km)
k = np.clip(np.searchsorted(tc, t), 1, len(tc) - 1); k = np.where(np.abs(tc[k - 1] - t) < np.abs(tc[k] - t), k - 1, k)
okt = np.abs(tc[k] - t) <= 1800
jj = np.abs(cla[None, :] - la[:, None]).argmin(1); ii = np.abs(clo[None, :] - lo[:, None]).argmin(1)
dkm = np.hypot((cla[jj] - la) * 111.19, (clo[ii] - lo) * 111.19 * np.cos(np.radians(38.1)))
u = cu[k, jj, ii]; v = cv[k, jj, ii]
okc = okt & (dkm < 3) & np.isfinite(u)
a = np.radians(dirc); wx, wy = -np.sin(a), -np.cos(a)            # verso dove va il vento
cw = u * wx + v * wy                                               # corrente lungo il vento (m/s)
reg = np.where((D0 >= 315) | (D0 < 45), 0, np.where((D0 >= 135) & (D0 < 225), 1, 2))
ok = okc & np.isfinite(vero) & (vero > 1) & (mod > 3)
y = np.log(vero / mod)
# media per cella e regime
key = (j * 100 + i) * 3 + reg
u_, inv = np.unique(key[ok], return_inverse=True)
m_ = np.bincount(inv, y[ok]) / np.bincount(inv)
res = np.full(len(y), np.nan); res[np.where(ok)[0]] = y[ok] - m_[inv]
x = -cw / mod
print(f"righe con corrente: {ok.sum()}, passaggi {len(np.unique(sc[ok]))}; |corrente| mediana {np.median(np.hypot(u,v)[ok])*1.944:.2f} nodi, 95% {np.percentile(np.hypot(u,v)[ok],95)*1.944:.2f}")
def retta(m, nome):
    xx, yy = x[m], res[m]
    b = np.sum((xx - xx.mean()) * (yy - yy.mean())) / np.sum((xx - xx.mean()) ** 2)
    r = np.corrcoef(xx, yy)[0, 1]
    # contro / con la corrente, corrente oltre 1 nodo
    forte = np.abs(cw[m]) > 0.5144
    contro = res[m][forte & (cw[m] < 0)].mean() if (forte & (cw[m] < 0)).sum() > 30 else np.nan
    con = res[m][forte & (cw[m] > 0)].mean() if (forte & (cw[m] > 0)).sum() > 30 else np.nan
    print(f"  {nome:28s} righe {m.sum():6d}  pendenza {b:+.2f} (1 = il radar vede il vento sull'acqua)  correlazione {r:+.2f}   "
          f"corrente oltre 1 nodo: contro il vento {100*(np.exp(contro)-1):+.1f}%, a favore {100*(np.exp(con)-1):+.1f}%")
for nome, m in (("tutto", ok), ("vento da nord", ok & (reg == 0)), ("vento da sud", ok & (reg == 1)), ("entro 3 km", ok & (dist < 3)),
                ("3-10 km", ok & (dist >= 3) & (dist < 10)), ("stretta (oltre 38,15)", ok & (la > 38.15) & (lo > 15.55)),
                ("vento oltre 6 m/s", ok & (mod > 6)), ("vento 3-6 m/s", ok & (mod <= 6))):
    if m.sum() > 200:
        retta(m, nome)
