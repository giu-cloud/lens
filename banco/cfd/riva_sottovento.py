"""Il radar sottovento ai rilievi (vento da terra, rilievo alto entro 3 km) a 1,5-3 km dalla riva: Lens 0.54 contro la
risposta piena dell'atlante, area per area (tabella dello studio statistico). Percorsi in TABELLA, VERO, CAMPI, ATLANTI.
"""
import os
import numpy as np, sys, json
sys.path.insert(0,'banco/cfd'); sys.path.insert(0,'banco/strumenti')
from riva_radar import risposta_su_radar, peso, AREE14
T=np.load(os.environ.get('TABELLA','/home/claude/lavoro/tabella6.npz'), mmap_mode='r'); V=np.load(os.environ.get('VERO','/home/claude/lavoro/radar_vero_tabella6.npy'), mmap_mode='r')
S='/home/claude/dati/shift'
area_all=np.asarray(T['area'])
tot={}
for nome in ["palermo","messina","napoli","genova","cagliari","egadi","eolie","malta","toscano","maddalena-bonifacio"]:
    m=area_all==AREE14.index(nome)
    C=np.load(f'/home/claude/lavoro/base3/campi_{nome}.npz')
    R=risposta_su_radar(f'{S}/atlante_{nome}.lna',C['lat'].astype(float),C['lon'].astype(float))
    j,i=T['j'][m].astype(int),T['i'][m].astype(int)
    vero=V[m]; mod=T['mod'][m]; q=T['q_lens054'][m]; D0=T['D0'][m]; dist=T['dist'][m]; terra=T['terra'][m]; ang=T['ang'][m]; S0=T['S0'][m]
    x=(D0%360)/22.5; k0=np.floor(x).astype(int)%16; w=x-np.floor(x); k1=(k0+1)%16
    r=(1-w)*R[k0,j,i]+w*R[k1,j,i]
    ok=np.isfinite(r)&np.isfinite(vero)&(mod>2)&(q>0)&(r>0)&(S0>4)
    base=mod*q
    for lab,cond in (("sottovento (terra a monte <3 km, rilievo >4°)",(terra<3)&(ang>4)),("sopravvento/lungo costa",(terra>=10))):
        for fascia,(lo,hi) in (("1.5-2",(1.5,2)),("2-3",(2,3))):
            f=ok&cond&(dist>=lo)&(dist<hi)
            if f.sum()<100: continue
            def mis(p):
                e=(p[f]-vero[f])/0.5144; return e.mean()
            key=(lab,fascia); t=tot.setdefault(key,{"n":0,"lens":0,"pieno":0,"riva05":0,"radar_su_lens":0})
            n=f.sum(); t["n"]+=n; t["lens"]+=mis(base)*n; t["pieno"]+=mis(base*r)*n; t["riva05"]+=mis(base*r**peso(dist,0.5,2))*n
            t["radar_su_lens"]+=np.median(vero[f]/base[f])*n
            print(f"{nome:20s} {lab[:12]} {fascia}: n {n:6d} radar/lens {np.median(vero[f]/base[f]):.2f} atlante {np.median(r[f]):.2f} | errore medio lens {mis(base):+.2f} pieno {mis(base*r):+.2f} kn")
for k,t in tot.items():
    n=t["n"]; print(k, n, {kk:round(v/n,3) for kk,v in t.items() if kk!="n"})
