"""Le "schermate" dello Stretto (26/11): Modelli | A | B | C prima (30 km) | C ora (5 km) | radar vero, per scene vere della
cassaforte, con i colori del vento dell'app (MapColors.wind). Campi da Schermate.kt.
Uso: python3 figura_schermate.py <scene.bin> <scene,...> <cartella png>"""
import sys, datetime as dt, numpy as np, matplotlib
from zoneinfo import ZoneInfo
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

STOP = [(0, '#7FB2E5'), (5, '#4CAF50'), (10, '#C0CA33'), (15, '#FFB300'), (20, '#FF7043'), (25, '#E53935'), (32, '#8E24AA')]
VENTO = LinearSegmentedColormap.from_list('app', [(k / 32, c) for k, c in STOP])
TERRA, COSTA, MARE = '#F3F0E8', '#8EA2B2', '#CFE3EC'
DIFF = LinearSegmentedColormap.from_list('diff', ['#2B6CB0', '#EDEDED', '#C0392B'])
fm = lambda x: f"{x:+.1f}".replace('.', ',')
KN = 1.943844


def leggi(p):
    b = open(p, 'rb').read()
    assert b[:4] == b'SCH1'
    ns, nx, ny = np.frombuffer(b, '<i4', 3, 4)
    box = np.frombuffer(b, '<f8', 4, 16)
    a = np.frombuffer(b, '<f4', offset=48).reshape(ns, 11, ny, nx)
    return box, a


def radar(scena):
    T = np.load('/home/claude/lavoro/tabella6.npz', mmap_mode='r'); sel = np.asarray(T['area']) == 3
    g = lambda k: np.asarray(T[k][sel])
    vero = np.load('/home/claude/lavoro/radar_vero_tabella6.npy', mmap_mode='r')[sel]
    C = np.load('/home/claude/lavoro/v_scia/campi_messina.npz'); rlat, rlon = C['lat'], C['lon']
    m = g('scena') == scena
    R = np.full((len(rlat), len(rlon)), np.nan); D = np.full_like(R, np.nan)
    R[g('j')[m].astype(int), g('i')[m].astype(int)] = vero[m]
    D[g('j')[m].astype(int), g('i')[m].astype(int)] = g('dist')[m]
    return rlat, rlon, R, D, int(g('t')[m][0])


def bordi(c):
    d = np.diff(c); return np.r_[c[0] - d[0] / 2, c[:-1] + d / 2, c[-1] + d[-1] / 2]


def main():
    box, A = leggi(sys.argv[1]); scene = [int(s) for s in sys.argv[2].split(',')]; cart = sys.argv[3]
    S, N, O, E = box; ny, nx = A.shape[2:]
    lat = N - (np.arange(ny) + 0.5) * (N - S) / ny; lon = O + (np.arange(nx) + 0.5) * (E - O) / nx
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    for k, sc in enumerate(scene):
        f = A[k]; mare = f[8] > 0.5
        sA = np.hypot(f[3], f[4])
        pann = [("Modelli (ECMWF)", f[0], f[1], f[2]), ("A, rilievo pieno", sA, f[3], f[4]),
                ("B, radar cella per cella", sA * f[5], f[3], f[4]), ("C prima, livello su 30 km", sA * f[6], f[3], f[4]),
                ("C ora, livello su 5 km", sA * f[7], f[3], f[4])]
        rlat, rlon, R, Dst, t = radar(sc)
        # confronto con il radar: media dei punti a 100 m dentro ogni cella del radar
        dj = abs(rlat[1] - rlat[0]); di = abs(rlon[1] - rlon[0])
        J = np.rint((LA - rlat[0]) / (rlat[1] - rlat[0])).astype(int); I = np.rint((LO - rlon[0]) / (rlon[1] - rlon[0])).astype(int)
        ok = (J >= 0) & (J < len(rlat)) & (I >= 0) & (I < len(rlon)) & mare
        def su_radar(s):
            q = ok & np.isfinite(s); key = J[q] * len(rlon) + I[q]
            sm = np.bincount(key, s[q], len(rlat) * len(rlon)); n = np.bincount(key, None, len(rlat) * len(rlon))
            out = np.where(n > 20, sm / np.maximum(n, 1), np.nan).reshape(len(rlat), len(rlon)); return out
        d = dt.datetime.fromtimestamp(t, dt.timezone.utc); loc = d.astimezone(ZoneInfo('Europe/Rome'))
        run = d.replace(hour=0 if d.hour < 12 else 12, minute=0)
        fig, AX = plt.subplots(2, 6, figsize=(30, 15.2), gridspec_kw={"hspace": 0.22}); ax = AX[0]
        for a in AX.ravel():
            a.set_facecolor(TERRA); a.set_xlim(O, E); a.set_ylim(S, N); a.set_aspect(1 / np.cos(np.radians((S + N) / 2)))
            a.set_xticks([]); a.set_yticks([])
        passo = 15
        for a, (nome, s, u, v) in zip(ax[:5], pann):
            im = a.pcolormesh(lon, lat, np.where(mare, s * KN, np.nan), cmap=VENTO, vmin=0, vmax=32, shading='auto', rasterized=True)
            a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
            sel = np.zeros_like(mare); sel[passo // 2::passo, passo // 2::passo] = True; sel &= mare & np.isfinite(u)
            sp = np.hypot(u, v) + 1e-9
            a.quiver(LO[sel], LA[sel], (u / sp)[sel], (v / sp)[sel], scale=38, width=0.004, color='#1F2D3A', alpha=0.65)
            x = su_radar(s); q = np.isfinite(x) & np.isfinite(R)
            e = (x - R)[q] * KN; e3 = (x - R)[q & (Dst < 3)] * KN
            a.set_title(f"{nome}\nrispetto al radar: {fm(e.mean())} nodi in media\n(errore quadratico {fm(np.sqrt((e**2).mean()))[1:]})\nentro 3 km dalla costa: {fm(e3.mean())}", fontsize=11)
        mod_r = su_radar(f[10])
        for k2, (nome, s_, u_, v_) in enumerate(pann):
            a = AX[1, k2]
            if k2 == 0:
                a.text(0.5, 0.5, "seconda riga:\nla differenza dai modelli\nlisciati su 4 km, come\nli riceve Lens\n(rosso = più vento,\nblu = meno)", transform=a.transAxes, ha='center', va='center', fontsize=13, color='#1F2D3A')
                a.set_facecolor('white'); a.axis('off'); continue
            imd = a.pcolormesh(lon, lat, np.where(mare, (s_ - f[10]) * KN, np.nan), cmap=DIFF, vmin=-10, vmax=10, shading='auto', rasterized=True)
            a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
            a.set_title(f"{nome.split(',')[0]} meno i modelli", fontsize=11)
        a = AX[1, 5]
        a.pcolormesh(lon, lat, np.where(mare, 1.0, np.nan), cmap=ListedColormap([MARE]), shading='auto', rasterized=True)
        a.pcolormesh(bordi(rlon), bordi(rlat), (R - mod_r) * KN, cmap=DIFF, vmin=-10, vmax=10, shading='flat', rasterized=True)
        a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
        a.set_title("Radar meno i modelli", fontsize=11)
        cb2 = fig.colorbar(imd, ax=AX[1].tolist(), shrink=0.75, pad=0.01, ticks=[-10, -5, 0, 5, 10]); cb2.set_label('nodi in più o in meno')
        a = ax[5]
        a.pcolormesh(lon, lat, np.where(mare, 1.0, np.nan), cmap=ListedColormap([MARE]), shading='auto', rasterized=True)
        a.pcolormesh(bordi(rlon), bordi(rlat), R * KN, cmap=VENTO, vmin=0, vmax=32, shading='flat', rasterized=True)
        a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
        a.set_title(f"Radar vero (Sentinel-1)\ncelle di 1 km; il primo km dalla\ncosta non lo vede (azzurro)\nmedia {np.nanmean(R) * KN:.1f} nodi".replace('.', ','), fontsize=11)
        for a in [x for x in AX.ravel() if x.axison]:
            for nome, la_, lo_ in (("Messina", 38.19, 15.555), ("Reggio", 38.11, 15.65), ("Villa S.G.", 38.22, 15.637), ("Ganzirri", 38.26, 15.61)):
                a.plot(lo_, la_, 'o', ms=2.5, color='#1F2D3A')
            a.text(15.53, 38.19, "Messina", fontsize=7, ha='right', color='#1F2D3A'); a.text(15.67, 38.11, "Reggio", fontsize=7, color='#1F2D3A')
        cb = fig.colorbar(im, ax=ax.tolist(), shrink=0.75, pad=0.01, ticks=[0, 5, 10, 15, 20, 25, 30]); cb.set_label('nodi')
        largo = pann[0][1][mare & (LA > 38.28)]
        fig.suptitle(f"Stretto di Messina, {loc:%d/%m/%Y %H:%M} ora italiana ({d:%H:%M} UTC), passaggio del radar. Previsione ECMWF del {run:%d/%m %H} UTC. "
                     f"Settimana della cassaforte: la libreria del radar non l'ha mai vista.", fontsize=15, y=0.97)
        fig.savefig(f"{cart}/stretto_schermata_{d:%Y%m%d_%H%M}.png", dpi=70, bbox_inches='tight')
        plt.close(fig)
        print(f"{cart}/stretto_schermata_{d:%Y%m%d_%H%M}.png")


main()
