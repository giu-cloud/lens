"""Strato marino a gravità ridotta sul rilievo: lo strato d'aria fresca sotto l'inversione che scorre attorno alla costa.

Equazioni delle acque basse per uno strato di spessore h sotto un'inversione (g' = g Δθ/θ), su griglia sfalsata (C):
  du/dt + u·∇u + f k×u = −g' ∇(h + b) − Cd |u| u / H + F   (attrito uguale su terra e mare: la scabrezza la fa
  il resto di Lens; qui contano solo pressione e inerzia, cioè il rilievo)
  dh/dt + ∇·(h u) = 0
b è il rilievo; dove supera 0,9 H fa da muro. F è la spinta della pressione a grande scala, scelta perché lontano dalla
costa il vento resti quello imposto (U∞, V∞) con h = H. Ai bordi una fascia di richiamo verso il vento imposto.
Il numero di Froude U/√(g'H) decide il regime: sotto 1 l'aria gira attorno ai rilievi (getti ai capi, ridossi), vicino
o sopra 1 accelera sui passi e scavalca le creste basse (getti nelle valli e negli stretti, salti idraulici).
Si integra fino a regime e si fa la media dell'ultima parte (le scie dietro le isole possono oscillare).

  from strato_marino import griglia, risolvi
  g = griglia(rilievo, sud, nord, ovest, est, passo_km=1.0)
  S, h = risolvi(g, U=0, V=-8, H=500, Fr=0.8)        # vento da nord a 8 m/s
"""
import math

import numba as nb
import numpy as np


def griglia(ril, s, n, w, e, passo_km=1.0, margine=0.35, liscio_km=8.0):
    """Griglia di calcolo con margine attorno al riquadro; rilievo medio sulle maglie (dal livello a 464 m)."""
    s0, n0, w0, e0 = s - margine, n + margine, w - margine * 1.25, e + margine * 1.25
    rl, ro, H = ril.griglia(0, s0, n0, w0, e0)
    H = np.maximum(H, 0).astype("f8")
    c = rl[1] - rl[0]
    latm = 0.5 * (s0 + n0)
    kmy = 111.19; kmx = 111.19 * math.cos(math.radians(latm))
    fy = max(1, int(round(passo_km / (c * kmy)))); fx = max(1, int(round(passo_km / (c * kmx))))
    ny, nx = H.shape[0] // fy, H.shape[1] // fx
    Hb = H[:ny * fy, :nx * fx].reshape(ny, fy, nx, fx)
    terra = (Hb > 0).mean(axis=(1, 3))
    # maglie di mare a livello del mare; per quelle di terra la quota media della sola parte di terra
    b = np.where(terra >= 0.5, np.where(Hb > 0, Hb, 0).sum(axis=(1, 3)) / np.maximum((Hb > 0).sum(axis=(1, 3)), 1), 0.0)
    lat = rl[0] + (np.arange(ny) * fy + (fy - 1) / 2) * c
    lon = ro[0] + (np.arange(nx) * fx + (fx - 1) / 2) * c
    dx = fx * c * kmx * 1000; dy = fy * c * kmy * 1000
    # il rilievo come lo vede il modello globale (ECMWF a 9 km, risoluzione efficace di qualche decina di km)
    from scipy.ndimage import gaussian_filter
    hm = (Hb.mean(axis=(1, 3)))
    b_l = gaussian_filter(hm, (liscio_km * 1000 / dy, liscio_km * 1000 / dx), mode="nearest")
    t_l = gaussian_filter(terra, (liscio_km * 1000 / dy, liscio_km * 1000 / dx), mode="nearest")
    return dict(lat=lat, lon=lon, b=b, terra=terra, dx=dx, dy=dy, f=2 * 7.292e-5 * math.sin(math.radians(latm)),
                b_liscio=np.where(t_l >= 0.5, b_l / np.maximum(t_l, 1e-3), 0.0), terra_liscia=t_l)


@nb.njit(cache=True)
def _d3(f0, fm1, fm2, fp1, fp2, c, d):
    """Derivata con differenze sopravvento del terz'ordine."""
    if c >= 0:
        return (2 * fp1 + 3 * f0 - 6 * fm1 + fm2) / (6 * d)
    return (-fp2 + 6 * fp1 - 3 * f0 - 2 * fm1) / (6 * d)


@nb.njit(cache=True)
def _tendenze(h, u, v, b, muro, cd, rich, U, V, Hm, gp, f, dx, dy, Fx, Fy, nu, du, dv, dh):
    ny, nx = h.shape
    HMIN = 10.0
    # u sulle facce x: (ny, nx+1); faccia i tra le celle i-1 e i
    for j in range(ny):
        for i in range(nx + 1):
            du[j, i] = 0.0
            if i == 0 or i == nx:
                continue
            if muro[j, i - 1] or muro[j, i]:
                continue
            uu = u[j, i]
            vv = 0.25 * (v[j, i - 1] + v[j, i] + v[j + 1, i - 1] + v[j + 1, i])
            # avvezione
            if 2 <= i <= nx - 2:
                ax = _d3(u[j, i], u[j, i - 1], u[j, i - 2], u[j, i + 1], u[j, i + 2], uu, dx)
            else:
                ax = (u[j, i] - u[j, i - 1]) / dx if uu >= 0 else (u[j, i + 1] - u[j, i]) / dx
            if 2 <= j <= ny - 3:
                ay = _d3(u[j, i], u[j - 1, i], u[j - 2, i], u[j + 1, i], u[j + 2, i], vv, dy)
            elif j == 0:
                ay = (u[j + 1, i] - u[j, i]) / dy
            elif j == ny - 1:
                ay = (u[j, i] - u[j - 1, i]) / dy
            else:
                ay = (u[j, i] - u[j - 1, i]) / dy if vv >= 0 else (u[j + 1, i] - u[j, i]) / dy
            eta1 = h[j, i] + b[j, i]; eta0 = h[j, i - 1] + b[j, i - 1]
            hf = Hm
            sp = math.sqrt(uu * uu + vv * vv)
            c = 0.5 * (cd[j, i] + cd[j, i - 1])
            lap = 0.0
            if 1 <= j <= ny - 2:
                lap = (u[j, i + 1] - 2 * uu + u[j, i - 1]) / (dx * dx) + (u[j + 1, i] - 2 * uu + u[j - 1, i]) / (dy * dy)
            r = 0.5 * (rich[j, i] + rich[j, i - 1])
            du[j, i] = (-uu * ax - vv * ay + f * vv - gp * (eta1 - eta0) / dx - c * sp * uu / hf + Fx + nu * lap
                        - r * (uu - U))
    # v sulle facce y: (ny+1, nx); faccia j tra le celle j-1 e j
    for j in range(ny + 1):
        for i in range(nx):
            dv[j, i] = 0.0
            if j == 0 or j == ny:
                continue
            if muro[j - 1, i] or muro[j, i]:
                continue
            vv = v[j, i]
            uu = 0.25 * (u[j - 1, i] + u[j - 1, i + 1] + u[j, i] + u[j, i + 1])
            if 2 <= j <= ny - 2:
                ay = _d3(v[j, i], v[j - 1, i], v[j - 2, i], v[j + 1, i], v[j + 2, i], vv, dy)
            else:
                ay = (v[j, i] - v[j - 1, i]) / dy if vv >= 0 else (v[j + 1, i] - v[j, i]) / dy
            if 2 <= i <= nx - 3:
                ax = _d3(v[j, i], v[j, i - 1], v[j, i - 2], v[j, i + 1], v[j, i + 2], uu, dx)
            elif i == 0:
                ax = (v[j, i + 1] - v[j, i]) / dx
            elif i == nx - 1:
                ax = (v[j, i] - v[j, i - 1]) / dx
            else:
                ax = (v[j, i] - v[j, i - 1]) / dx if uu >= 0 else (v[j, i + 1] - v[j, i]) / dx
            eta1 = h[j, i] + b[j, i]; eta0 = h[j - 1, i] + b[j - 1, i]
            hf = Hm
            sp = math.sqrt(uu * uu + vv * vv)
            c = 0.5 * (cd[j, i] + cd[j - 1, i])
            lap = 0.0
            if 1 <= i <= nx - 2:
                lap = (v[j + 1, i] - 2 * vv + v[j - 1, i]) / (dy * dy) + (v[j, i + 1] - 2 * vv + v[j, i - 1]) / (dx * dx)
            r = 0.5 * (rich[j, i] + rich[j - 1, i])
            dv[j, i] = (-uu * ax - vv * ay - f * uu - gp * (eta1 - eta0) / dy - c * sp * vv / hf + Fy + nu * lap
                        - r * (vv - V))
    # continuità, flussi sopravvento
    for j in range(ny):
        for i in range(nx):
            dh[j, i] = 0.0
            if muro[j, i]:
                continue
            fl = 0.0
            # facce x
            ue = u[j, i + 1]; uw = u[j, i]
            he = h[j, i] if ue >= 0 else (h[j, i + 1] if i + 1 < nx else h[j, i])
            hw = (h[j, i - 1] if i - 1 >= 0 else h[j, i]) if uw >= 0 else h[j, i]
            fl -= (ue * he - uw * hw) / dx
            vn = v[j + 1, i]; vs = v[j, i]
            hn = h[j, i] if vn >= 0 else (h[j + 1, i] if j + 1 < ny else h[j, i])
            hs = (h[j - 1, i] if j - 1 >= 0 else h[j, i]) if vs >= 0 else h[j, i]
            fl -= (vn * hn - vs * hs) / dy
            href = max(Hm - b[j, i], 0.0)
            dh[j, i] = fl - rich[j, i] * (h[j, i] - href)


@nb.njit(cache=True)
def _integra(h, u, v, b, muro, cd, rich, U, V, Hm, gp, f, dx, dy, Fx, Fy, nu, dt, nsteps, nmedia):
    ny, nx = h.shape
    du = np.zeros_like(u); dv = np.zeros_like(v); dh = np.zeros_like(h)
    u1 = u.copy(); v1 = v.copy(); h1 = h.copy()
    Su = np.zeros((ny, nx)); Sv = np.zeros((ny, nx)); Ss = np.zeros((ny, nx)); Sh = np.zeros((ny, nx))
    nm = 0
    for k in range(nsteps):
        # Runge-Kutta a tre passi (Wicker e Skamarock)
        for st in range(3):
            a = dt / 3.0 if st == 0 else (dt / 2.0 if st == 1 else dt)
            _tendenze(h1, u1, v1, b, muro, cd, rich, U, V, Hm, gp, f, dx, dy, Fx, Fy, nu, du, dv, dh)
            for j in range(ny):
                for i in range(nx + 1):
                    u1[j, i] = u[j, i] + a * du[j, i]
            for j in range(ny + 1):
                for i in range(nx):
                    v1[j, i] = v[j, i] + a * dv[j, i]
            for j in range(ny):
                for i in range(nx):
                    h1[j, i] = max(h[j, i] + a * dh[j, i], 1.0)
        u[:, :] = u1; v[:, :] = v1; h[:, :] = h1
        if k >= nsteps - nmedia:
            for j in range(ny):
                for i in range(nx):
                    uc = 0.5 * (u[j, i] + u[j, i + 1]); vc = 0.5 * (v[j, i] + v[j + 1, i])
                    Su[j, i] += uc; Sv[j, i] += vc; Ss[j, i] += math.sqrt(uc * uc + vc * vc); Sh[j, i] += h[j, i]
            nm += 1
    return Su / nm, Sv / nm, Ss / nm, Sh / nm


def risolvi(g, U, V, H=500.0, Fr=0.8, cd_mare=1.2e-3, cd_terra=1.2e-3, quota_muro=0.9, nu=50.0, attraversamenti=2.5,
            fascia=14, tau=600.0, liscio=False):
    """Vento dello strato (u, v, velocità media, spessore) sulla griglia g, con vento imposto (U, V) m/s.
    liscio=True: sul rilievo come lo vede il modello globale."""
    b = g["b_liscio"] if liscio else g["b"]; ny, nx = b.shape
    terra_g = g["terra_liscia"] if liscio else g["terra"]
    W = max(math.hypot(U, V), 0.5)
    gp = (W / Fr) ** 2 / H
    muro = b >= quota_muro * H
    cd = np.where(terra_g > 0.5, cd_terra, cd_mare)
    jj, ii = np.indices(b.shape)
    d = np.minimum(np.minimum(jj, ny - 1 - jj), np.minimum(ii, nx - 1 - ii)).astype("f8")
    rich = np.where(d < fascia, ((fascia - d) / fascia) ** 2 / tau, 0.0)
    f = g["f"]
    Fx = -f * V + cd_mare * W * U / H
    Fy = f * U + cd_mare * W * V / H
    h = np.where(muro, 1.0, np.maximum(H - b, 1.0))
    u = np.full((ny, nx + 1), float(U)); v = np.full((ny + 1, nx), float(V))
    for j in range(ny):
        for i in range(nx + 1):
            if (i > 0 and muro[j, i - 1]) or (i < nx and muro[j, i]):
                u[j, i] = 0.0
    for j in range(ny + 1):
        for i in range(nx):
            if (j > 0 and muro[j - 1, i]) or (j < ny and muro[j, i]):
                v[j, i] = 0.0
    c = math.sqrt(gp * H)
    dt = 0.45 * min(g["dx"], g["dy"]) / (2.0 * W + c)
    L = max(nx * g["dx"], ny * g["dy"])
    tempo = attraversamenti * L / W + 3 * 3600
    n = int(tempo / dt)
    nmedia = int(0.35 * n)
    Su, Sv, Ss, Sh = _integra(h, u, v, b.astype("f8"), muro, cd, rich, float(U), float(V), float(H), gp, f, g["dx"],
                              g["dy"], Fx, Fy, nu, dt, n, nmedia)
    Ss = np.where(muro, np.nan, Ss)
    return dict(u=Su, v=Sv, S=Ss, h=Sh, muro=muro, passi=n, dt=dt, gp=gp)


def lettura_mare(g, campo, lat_t, lon_t, passo=0.25, terra_nodo=None):
    """Il campo letto come lo legge Shift dal modello a 0,25°: solo i punti di mare della griglia, pesi bilineari
    rinormalizzati; se i quattro sono a terra, il punto di mare più vicino entro due maglie.
    campo: (ny, nx) sulla griglia g; lat_t, lon_t: vettori della griglia di arrivo. Torna (len(lat_t), len(lon_t))."""
    lat, lon = g["lat"], g["lon"]
    la0 = math.floor(lat.min() / passo) * passo; lo0 = math.floor(lon.min() / passo) * passo
    nla = int(math.ceil((lat.max() - la0) / passo)) + 1; nlo = int(math.ceil((lon.max() - lo0) / passo)) + 1
    NL = la0 + np.arange(nla) * passo; NO = lo0 + np.arange(nlo) * passo
    jn = np.clip(np.round((NL - lat[0]) / (lat[1] - lat[0])).astype(int), 0, len(lat) - 1)
    in_ = np.clip(np.round((NO - lon[0]) / (lon[1] - lon[0])).astype(int), 0, len(lon) - 1)
    val = campo[np.ix_(jn, in_)]
    if terra_nodo is None:
        terra_nodo = g["terra"][np.ix_(jn, in_)] > 0.5
    terra_nodo = terra_nodo | ~np.isfinite(val)
    out = np.full((len(lat_t), len(lon_t)), np.nan)
    for a, la in enumerate(lat_t):
        fj = (la - la0) / passo; j0 = int(math.floor(fj)); y = fj - j0
        for o, lo in enumerate(lon_t):
            fi = (lo - lo0) / passo; i0 = int(math.floor(fi)); x = fi - i0
            pts = [(j0, i0, (1 - x) * (1 - y)), (j0, i0 + 1, x * (1 - y)), (j0 + 1, i0, (1 - x) * y), (j0 + 1, i0 + 1, x * y)]
            s = 0.0; ws = 0.0
            for j, i, w in pts:
                if 0 <= j < nla and 0 <= i < nlo and not terra_nodo[j, i] and w > 0:
                    s += w * val[j, i]; ws += w
            if ws > 1e-6:
                out[a, o] = s / ws
                continue
            best = None; bd = 1e9
            for j in range(j0 - 2, j0 + 4):
                for i in range(i0 - 2, i0 + 4):
                    if 0 <= j < nla and 0 <= i < nlo and not terra_nodo[j, i]:
                        dd = (j - fj) ** 2 + (i - fi) ** 2
                        if dd < bd:
                            bd = dd; best = val[j, i]
            if best is not None:
                out[a, o] = best
    return out


def al_punto(g, campo, lat_t, lon_t):
    """Il campo della griglia g nelle celle (lat_t × lon_t), bilineare."""
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((g["lat"], g["lon"]), campo, bounds_error=False, fill_value=np.nan)
    LA, LO = np.meshgrid(lat_t, lon_t, indexing="ij")
    return f((LA, LO))
