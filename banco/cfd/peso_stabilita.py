"""Il peso della CFD secondo forza del vento al largo e stabilità (dT925), sulle righe di taratura_radar.py --npz.

  python3 peso_stabilita.py <righe .npz> <area>
"""
import os
import numpy as np, sys
f=sys.argv[1]; area=sys.argv[2]
AREE14=["egadi","eolie","palermo","messina","malta","napoli","cagliari","maddalena-bonifacio","toscano","genova","pontine","corsica-nord","monaco","cicladi"]
z=np.load(f); T=np.load(os.environ.get('TABELLA','/home/claude/lavoro/tabella6.npz'), mmap_mode='r')
m=T['area']==AREE14.index(area)
dT925=T['dT925'][m]; v925=T['v925'][m]
ok=z['ok']; y=z['y']; lL=z['lL']; lC=z['lC']; vero=z['vero']; mod=z['mod']; S0=z['S0']; cass=z['cass']; j=z['j']; i=z['i']; dist=z['dist']
ok=ok&np.isfinite(dT925)
cel=j*1000+i
def sist(e,mask,cel):
    ids,inv=np.unique(cel[mask],return_inverse=True); n=np.bincount(inv); me=np.bincount(inv,e[mask])/n
    return np.sqrt(np.sum(n*me**2)/n.sum())
tr=ok&~cass; te=ok&cass
# peso della CFD: funzione liscia di forza al largo e stabilità
ws=np.clip((S0-3)/6,0,1)                 # 0 sotto 3 m/s, 1 sopra 9 m/s
st=np.clip((4-dT925)/4,0,1)              # 1 se dT925<=0 (neutro/instabile), 0 se >=4 K (stabile)
def valuta(X,nome):
    X1=np.column_stack([np.ones(tr.sum())]+[x[tr] for x in X]); c=np.linalg.lstsq(X1,y[tr],rcond=None)[0]
    p=c[0]+sum(ci*x[te] for ci,x in zip(c[1:],X)); e=(mod[te]*np.exp(p)-vero[te])/0.5144; v=dist[te]<3
    print(f"{area} {nome:40s} {np.round(c,2).tolist()} rmse {np.sqrt(np.mean(e**2)):.3f} sist {sist(e,np.ones(te.sum(),bool),cel[te]):.3f} <3km sist {sist(e,v,cel[te]):.3f}")
valuta([],"modello"); valuta([lL],"lens"); valuta([lC],"cfd"); valuta([lL,lC],"lens+cfd")
valuta([lL,lC*ws],"lens+cfd*forza"); valuta([lL,lC*st],"lens+cfd*neutralità"); valuta([lL,lC*ws*st],"lens+cfd*forza*neutralità")
valuta([lC*ws*st],"cfd*forza*neutralità")
for nome,q in (("forza",ws),("neutralita",st)):
    print(area, nome, "media", q[ok].mean().round(2))
