"""Figure della revisione del 27/11 sera: Modelli | C oggi (livello su 30 km) | C proposto (livello su 10 km) | radar vero,
per passaggi veri del radar nelle settimane della cassaforte, colori del vento dell'app. Campi da SchermateRev.kt.
Regola di scelta dei passaggi, fissata prima di guardare: tra quelli della cassaforte con il modello ad almeno 8 nodi di
mediana, i due con più celle del radar dentro il riquadro dell'area fine.
Uso: python3 figura_rev.py <area> <cartella png>      (aree: palermo eolie corsica-nord dalmazia)
"""
import sys, os, subprocess, datetime as dt, numpy as np, matplotlib
from zoneinfo import ZoneInfo
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cfd"))
from riferimenti import AREE, L, D, MAPPE

STOP = [(0, '#7FB2E5'), (5, '#4CAF50'), (10, '#C0CA33'), (15, '#FFB300'), (20, '#FF7043'), (25, '#E53935'), (32, '#8E24AA')]
VENTO = LinearSegmentedColormap.from_list('app', [(k / 32, c) for k, c in STOP])
TERRA, COSTA, MARE = '#F3F0E8', '#8EA2B2', '#CFE3EC'
DIFF = LinearSegmentedColormap.from_list('diff', ['#2B6CB0', '#EDEDED', '#C0392B'])
fm = lambda x: f"{x:+.1f}".replace('.', ',')
KN = 1.943844
FINE = {"palermo": (f"{D}/atlanti_fini2/atlante_palermo-fine.lna", "Golfo di Palermo"), "eolie": (f"{D}/atlanti_fini2/atlante_eolie-centro.lna", "Eolie, Lipari Salina Vulcano"),
        "corsica-nord": (f"{D}/atlanti_fini3/atlante_corsica-capocorso.lna", "Capo Corso"), "dalmazia": (f"{D}/atlanti_fini3/atlante_dalmazia-split.lna", "Dalmazia, Spalato")}
CP = "/tmp/lens-motore-build/motore.jar:/tmp/rev/sch.jar:/opt/kotlin2021/kotlinc/lib/kotlin-stdlib.jar"
KA, KB = "30", "10"


def bordi(c):
    d = np.diff(c); return np.r_[c[0] - d[0] / 2, c[:-1] + d / 2, c[-1] + d[-1] / 2]


def main():
    nome, cart = sys.argv[1], sys.argv[2]; c = AREE[nome]; fine, titolo = FINE[nome]
    T = np.load(f"{L}/{c['tab']}.npz", mmap_mode="r"); sel = np.asarray(T["area"]) == c["idx"]
    g = lambda k: np.asarray(T[k][sel])
    vero = np.load(f"{L}/radar_vero_{c['tab']}.npy", mmap_mode="r")[sel]
    G = np.load(c["griglia"]); rlat, rlon = np.asarray(G["lat"]), np.asarray(G["lon"])
    sc, t, j, i, mod, dist, v925, d925 = g("scena"), g("t"), g("j").astype(int), g("i").astype(int), g("mod"), g("dist"), g("v925"), g("d925")
    la, lo = rlat[j], rlon[i]; b = c["box"]
    dentro = (la >= b[0]) & (la <= b[1]) & (lo >= b[2]) & (lo <= b[3]) & np.isfinite(vero) & ((t // (7 * 86400)) % 5 == 2)
    cand = []
    for s in np.unique(sc[dentro]):
        m = dentro & (sc == s); tt = int(t[m][0]); d = dt.datetime.fromtimestamp(tt, dt.timezone.utc)
        run = d.replace(hour=0 if d.hour < 12 else 12, minute=0, second=0); pw = f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2"
        if os.path.exists(pw) and np.median(mod[m]) * KN >= 8:
            sp = np.nanmedian(v925[m]); dd = np.radians(np.nanmedian(d925[m]))
            cand.append((m.sum(), s, tt, pw, -sp * np.sin(dd), -sp * np.cos(dd)))
    cand.sort(reverse=True); scelti = cand[:2]
    usc = f"/tmp/rev/sch_{nome}.bin"
    cmd = ["java", "-Xmx5g", "-cp", CP, "SchermateRevKt", c["atl"], fine, c["lib"], MAPPE, usc, KA, KB]
    for n, s, tt, pw, u, v in scelti:
        cmd += [pw, str(tt * 1000), f"{u:.3f}" if np.isfinite(u) else "nan", f"{v:.3f}" if np.isfinite(v) else "nan"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:]); raise SystemExit(1)
    raw = open(usc, 'rb').read(); ns, nx, ny = np.frombuffer(raw, '<i4', 3, 4); box = np.frombuffer(raw, '<f8', 4, 16)
    A = np.frombuffer(raw, '<f4', offset=48).reshape(ns, 11, ny, nx)
    S, N, O, E = box
    lat = N - (np.arange(ny) + 0.5) * (N - S) / ny; lon = O + (np.arange(nx) + 0.5) * (E - O) / nx
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    for k, (n, s, tt, pw, u9, v9) in enumerate(scelti):
        f = A[k]; mare = f[8] > 0.5; sA = np.hypot(f[3], f[4])
        pann = [("Modelli (ECMWF)", f[0], f[1], f[2]), ("C oggi: livello su 30 km", sA * f[6], f[3], f[4]),
                ("C proposto: livello su 10 km", sA * f[7], f[3], f[4])]
        m = sc == s
        R = np.full((len(rlat), len(rlon)), np.nan); Dst = np.full_like(R, np.nan)
        R[j[m], i[m]] = vero[m]; Dst[j[m], i[m]] = dist[m]
        J = np.rint((LA - rlat[0]) / (rlat[1] - rlat[0])).astype(int); I = np.rint((LO - rlon[0]) / (rlon[1] - rlon[0])).astype(int)
        ok = (J >= 0) & (J < len(rlat)) & (I >= 0) & (I < len(rlon)) & mare
        def su_radar(x):
            q = ok & np.isfinite(x); key = J[q] * len(rlon) + I[q]
            sm = np.bincount(key, x[q], len(rlat) * len(rlon)); nn = np.bincount(key, None, len(rlat) * len(rlon))
            return np.where(nn > 20, sm / np.maximum(nn, 1), np.nan).reshape(len(rlat), len(rlon))
        d = dt.datetime.fromtimestamp(tt, dt.timezone.utc); loc = d.astimezone(ZoneInfo('Europe/Rome'))
        run = d.replace(hour=0 if d.hour < 12 else 12, minute=0)
        asp = (N - S) / ((E - O) * np.cos(np.radians((S + N) / 2)))
        fig, AX = plt.subplots(2, 4, figsize=(26, 2 * 5.6 * asp + 3.2), gridspec_kw={"hspace": 0.2})
        for a in AX.ravel():
            a.set_facecolor(TERRA); a.set_xlim(O, E); a.set_ylim(S, N); a.set_aspect(1 / np.cos(np.radians((S + N) / 2))); a.set_xticks([]); a.set_yticks([])
        passo = max(8, nx // 28)
        mod_r = su_radar(f[10])
        for k2, (nm, s_, u, v) in enumerate(pann):
            a = AX[0, k2]
            im = a.pcolormesh(lon, lat, np.where(mare, s_ * KN, np.nan), cmap=VENTO, vmin=0, vmax=32, shading='auto', rasterized=True)
            a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
            q = np.zeros_like(mare); q[passo // 2::passo, passo // 2::passo] = True; q &= mare & np.isfinite(u)
            sp = np.hypot(u, v) + 1e-9
            a.quiver(LO[q], LA[q], (u / sp)[q], (v / sp)[q], scale=38, width=0.004, color='#1F2D3A', alpha=0.65)
            x = su_radar(s_); qq = np.isfinite(x) & np.isfinite(R)
            e = (x - R)[qq] * KN; e3 = (x - R)[qq & (Dst < 3)] * KN
            a.set_title(f"{nm}\nrispetto al radar: {fm(e.mean())} nodi in media\nerrore quadratico {fm(np.sqrt((e**2).mean()))[1:]}; entro 3 km dalla costa {fm(e3.mean()) if len(e3) else '-'}", fontsize=11)
            a = AX[1, k2]
            if k2 == 0:
                a.text(0.5, 0.5, "seconda riga:\nla differenza dai modelli\nlisciati su 4 km, come\nli riceve Lens\n(rosso = più vento,\nblu = meno)", transform=a.transAxes, ha='center', va='center', fontsize=13, color='#1F2D3A')
                a.set_facecolor('white'); a.axis('off'); continue
            imd = a.pcolormesh(lon, lat, np.where(mare, (s_ - f[10]) * KN, np.nan), cmap=DIFF, vmin=-10, vmax=10, shading='auto', rasterized=True)
            a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
            a.set_title(f"{nm.split(':')[0]} meno i modelli", fontsize=12)
        a = AX[0, 3]
        a.pcolormesh(lon, lat, np.where(mare, 1.0, np.nan), cmap=ListedColormap([MARE]), shading='auto', rasterized=True)
        a.pcolormesh(bordi(rlon), bordi(rlat), R * KN, cmap=VENTO, vmin=0, vmax=32, shading='flat', rasterized=True)
        a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
        a.set_title("Radar vero (Sentinel-1)\ncelle di 1 km; il primo chilometro e mezzo\ndalla costa non lo vede (azzurro)", fontsize=12)
        a = AX[1, 3]
        a.pcolormesh(lon, lat, np.where(mare, 1.0, np.nan), cmap=ListedColormap([MARE]), shading='auto', rasterized=True)
        a.pcolormesh(bordi(rlon), bordi(rlat), (R - mod_r) * KN, cmap=DIFF, vmin=-10, vmax=10, shading='flat', rasterized=True)
        a.contour(lon, lat, mare.astype(float), [0.5], colors=COSTA, linewidths=0.7)
        a.set_title("Radar meno i modelli", fontsize=12)
        fig.colorbar(im, ax=AX[0].tolist(), shrink=0.8, pad=0.01, ticks=[0, 5, 10, 15, 20, 25, 30]).set_label('nodi')
        fig.colorbar(imd, ax=AX[1].tolist(), shrink=0.8, pad=0.01, ticks=[-10, -5, 0, 5, 10]).set_label('nodi in più o in meno')
        fig.suptitle(f"{titolo}, {loc:%d/%m/%Y %H:%M} ora italiana ({d:%H:%M} UTC), passaggio del radar. Previsione ECMWF del {run:%d/%m %H} UTC. "
                     f"Settimana della cassaforte: la libreria del radar non l'ha mai vista.", fontsize=14, y=0.99)
        p = f"{cart}/rev_{nome}_{d:%Y%m%d_%H%M}.png"
        fig.savefig(p, dpi=62, bbox_inches='tight'); plt.close(fig); print(p)


main()
