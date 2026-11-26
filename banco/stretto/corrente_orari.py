"""Ore dei massimi di corrente di marea (montante = verso nord, scendente = verso sud) alla stretta, secondo il modello Copernicus MED.
Da confrontare con le tavole delle correnti dell'Istituto Idrografico della Marina. Ore italiane (UTC+1 d'inverno, +2 d'estate)."""
import numpy as np, datetime as dt, json, sys
Z = np.load('/home/claude/dati/correnti/correnti_stretto.npz')
tc, la, lo, v, u = Z['t'], Z['lat'], Z['lon'], Z['v'].astype(float), Z['u'].astype(float)
j, i = int(np.abs(la - 38.229).argmin()), int(np.abs(lo - 15.625).argmin())
s = v[:, j, i]                                           # verso nord, m/s (corrente totale: il "resto" è piccolo rispetto alla marea)
giorni = sys.argv[1:] or ["2026-11-22", "2026-11-23", "2026-11-24", "2026-11-25"]
out = {}
for g in giorni:
    d0 = dt.datetime.fromisoformat(g).replace(tzinfo=dt.timezone.utc).timestamp()
    m = (tc >= d0 - 7200) & (tc < d0 + 86400 + 7200); idx = np.where(m)[0]
    ev = []
    for k in idx[1:-1]:
        if (s[k] > s[k-1] and s[k] >= s[k+1] and s[k] > 0.15) or (s[k] < s[k-1] and s[k] <= s[k+1] and s[k] < -0.15):
            # vertice della parabola per l'ora fra i campioni
            a, b, c = s[k-1], s[k], s[k+1]; den = a - 2*b + c; f = 0.5 * (a - c) / den if den != 0 else 0.0
            t = tc[k] + f * 3600
            loc = dt.datetime.fromtimestamp(t, dt.timezone(dt.timedelta(hours=1)))
            if loc.date().isoformat() == g:
                ev.append((loc.strftime('%H:%M'), "montante" if b > 0 else "scendente", round(abs(b) * 1.944, 1)))
    out[g] = ev
    print(g, "  ".join(f"{h} {n} {kn} nodi" for h, n, kn in ev))
json.dump(out, open('/home/claude/lavoro/stretto/orari_corrente.json', 'w'), indent=1)
