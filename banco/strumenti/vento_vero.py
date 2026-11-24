"""Dal vento del radar (vento neutro equivalente a 10 m) al vento vero a 10 m, con la stabilità di ECMWF.

Il radar misura la scabrezza del mare, cioè lo sforzo del vento (u*), e la traduce nel vento che ci sarebbe a 10 m se
l'aria fosse neutra: U10N = u*/κ · ln(10/z0). Il vento vero è U10 = U10N − u*/κ · ψm(10/L). Con aria instabile (mare più
caldo dell'aria) ψm > 0 e il vento vero è più debole di quello del radar; con aria stabile ψm < 0 e il vero è più forte.
L, la lunghezza di Obukhov, viene dalla differenza di temperatura (virtuale) tra l'aria a 2 m e il mare: formule di
massa (Businger-Dyer, Paulson; Charnock per z0), qualche iterazione.

  from vento_vero import neutro_a_vero
  U10 = neutro_a_vero(U10N, dTam, sst, umid)
dTam = 2 m − mare (K), sst in °C, umid = 2 m − rugiada (K).
"""
import numpy as np

KAPPA = 0.4
G = 9.81


def _psi(zeta):
    zeta = np.asarray(zeta, dtype="f8")
    instab = zeta < 0
    x = np.power(np.maximum(1 - 16 * np.minimum(zeta, 0), 1.0), 0.25)
    psim_i = 2 * np.log((1 + x) / 2) + np.log((1 + x * x) / 2) - 2 * np.arctan(x) + np.pi / 2
    psih_i = 2 * np.log((1 + x * x) / 2)
    zs = np.clip(zeta, 0, 1.0)
    psim_s = -5 * zs; psih_s = -5 * zs
    return np.where(instab, psim_i, psim_s), np.where(instab, psih_i, psih_s)


def _qsat(Tc, p=1013.0):
    e = 6.112 * np.exp(17.67 * Tc / (Tc + 243.5))
    return 0.622 * e / (p - 0.378 * e)


def neutro_a_vero(un, dTam, sst, umid, iterazioni=6, ritorna_zeta=False):
    un = np.maximum(np.asarray(un, dtype="f8"), 0.5)
    sst = np.asarray(sst, dtype="f8"); dTam = np.asarray(dTam, dtype="f8"); umid = np.asarray(umid, dtype="f8")
    t2 = sst + dTam
    td = t2 - np.maximum(umid, 0)
    qa = _qsat(td); qs = 0.98 * _qsat(sst)
    th = t2 + 273.15
    dthv = (t2 + 0.0098 * 2 - sst) + 0.61 * th * (qa - qs)
    # u* e z0 dal vento neutro (Charnock), si risolve per punto fisso
    us = un * 0.035
    for _ in range(8):
        z0 = 0.011 * us * us / G + 0.11 * 1.5e-5 / np.maximum(us, 1e-3)
        us = KAPPA * un / np.log(10.0 / z0)
    z0t = 1.0e-4
    zeta10 = np.zeros_like(un)
    for _ in range(iterazioni):
        L = np.where(np.abs(zeta10) > 1e-9, 10.0 / np.where(np.abs(zeta10) > 1e-9, zeta10, 1.0), 1e9)
        _, psih2 = _psi(2.0 / L)
        tstar = KAPPA * dthv / (np.log(2.0 / z0t) - psih2)
        Lnew = th * us * us / (KAPPA * G * np.where(np.abs(tstar) > 1e-6, tstar, 1e-6))
        zeta10 = np.clip(10.0 / Lnew, -5, 1)
    psim10, _ = _psi(zeta10)
    u10 = np.maximum(un - us / KAPPA * psim10, 0.0)
    return (u10, zeta10) if ritorna_zeta else u10
