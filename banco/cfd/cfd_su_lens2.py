"""La CFD aggiunge qualcosa sopra Lens 2? Fattore della libreria senza cassaforte sulle righe della cassaforte, poi
Lens 2 da solo e Lens 2 con la CFD, stimati su metà delle settimane e provati sull'altra metà (due pieghe).

  python3 cfd_su_lens2.py <righe .npz> <area> <libreria senza cassaforte .lra>
"""
import os
import numpy as np, sys
sys.path.insert(0,'banco/strumenti')
from libreria_radar import Libreria
f, area, libp = sys.argv[1:4]
AREE14=["egadi","eolie","palermo","messina","malta","napoli","cagliari","maddalena-bonifacio","toscano","genova","pontine","corsica-nord","monaco","cicladi"]
T=np.load(os.environ.get('TABELLA','/home/claude/lavoro/tabella6.npz'), mmap_mode='r'); m=np.asarray(T['area'])==AREE14.index(area)
z=np.load(f)
L=Libreria(libp)
ok=z['ok']; vero=z['vero']; mod=z['mod']; lL=z['lL']; lC=z['lC']; dist=z['dist']; S0=z['S0']; D0=z['D0']; sc=z['scena']; cass=z['cass']; j=z['j']; i=z['i']
dirc=T['dir'][m]; v925=T['v925'][m]; d925=T['d925'][m]; ora=T['ora'][m]; dT925=T['dT925'][m]; t=T['t'][m]
cella=j*L.nx+i; pos={c:k for k,c in enumerate(L.celle)}
G=np.full(len(j),np.nan)
for s in np.unique(sc[ok&cass]):
    r=np.where((sc==s))[0]
    a=np.radians(D0[r[0]]); b=np.radians(d925[r[0]])
    # direzione del modello per cella della libreria: dalle righe della scena (le altre celle: D0)
    dc=np.full(L.nc, D0[r[0]]); idx=[pos.get(c,-1) for c in cella[r]]
    for k,ii in zip(r,idx):
        if ii>=0: dc[ii]=dirc[k]
    g=L.fattore(S0[r[0]]*np.sin(a),S0[r[0]]*np.cos(a),v925[r[0]]*np.sin(b),v925[r[0]]*np.cos(b),1.0 if ora[r[0]]>=12 else 0.0,dc)
    for k,ii in zip(r,idx):
        if ii>=0: G[k]=g[ii]
use=ok&cass&np.isfinite(G)&np.isfinite(dT925)
neu=np.clip((4-dT925)/4,0,1)
sett=(t//(7*86400))
cel=j*1000+i
def sist(e,mask):
    ids,inv=np.unique(cel[mask],return_inverse=True); n=np.bincount(inv); me=np.bincount(inv,e[mask])/n
    return np.sqrt(np.sum(n*me**2)/n.sum())
y=np.log(vero/mod)
out={}
for nome,X in (("lens2",[np.log(G)]),("lens2+cfd",[np.log(G),lC]),("lens2+cfd*neutralità",[np.log(G),lC*neu]),("lens2+lens054",[np.log(G),lL])):
    errs=np.zeros(use.sum()); prev=np.zeros(len(y))
    for fold in (0,1):
        trn=use&((sett//5)%2==fold); tst=use&((sett//5)%2!=fold)
        X1=np.column_stack([np.ones(trn.sum())]+[x[trn] for x in X]); c=np.linalg.lstsq(X1,y[trn],rcond=None)[0]
        prev[tst]=c[0]+sum(ci*x[tst] for ci,x in zip(c[1:],X))
        if fold==0: c0=c
    e=(mod*np.exp(prev)-vero)/0.5144
    print(f"{area} {nome:24s} coeff(piega 0) {np.round(c0,2).tolist()} rmse {np.sqrt(np.mean(e[use]**2)):.3f} sist {sist(e,use):.3f} <3km sist {sist(e,use&(dist<3)):.3f} (n {use.sum()})")
e=(mod-vero)/0.5144; print(f"{area} modello grezzo rmse {np.sqrt(np.mean(e[use]**2)):.3f} sist {sist(e,use):.3f}")
