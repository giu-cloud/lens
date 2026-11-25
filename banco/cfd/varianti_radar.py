"""Varianti della CFD contro il radar: media pesata sulle direzioni (incertezza della direzione del modello) e
lisciatura nello spazio; stima fuori dalla cassaforte, prova nella cassaforte.

  python3 varianti_radar.py <righe .npz di taratura_radar.py>
"""
import os
import numpy as np, sys
from scipy.ndimage import gaussian_filter
f=sys.argv[1]
z=np.load(f)
ok=z['ok']; y=z['y']; lL=z['lL']; vero=z['vero']; mod=z['mod']; D0=z['D0']; S0=z['S0']; sc=z['scena']; j=z['j']; i=z['i']; dist=z['dist']
Sr=z['Sr']; Soff=z['Soff']; da=z['da']
T=np.load(os.environ.get('TABELLA','/home/claude/lavoro/tabella6.npz'), mmap_mode='r')
cass=z['cass'] if 'cass' in z.files else None
def lc_var(sig_d, sig_s):
    R=Sr/Soff[:,None,None]
    if sig_s>0:
        Rs=np.empty_like(R)
        for k in range(len(da)):
            a=R[k]; m=np.isfinite(a); a0=np.where(m,a,0)
            Rs[k]=gaussian_filter(a0,sig_s)/np.maximum(gaussian_filter(m.astype(float),sig_s),1e-6); Rs[k][~m]=np.nan
        R=Rs
    d=((D0[:,None]-da[None,:]+540)%360-180)
    w=np.clip(1-np.abs(d)/22.5,0,None) if sig_d==0 else np.exp(-0.5*(d/sig_d)**2)
    Rr=R[:,j,i].T
    w=np.where(np.isfinite(Rr),w,0); w=w/np.maximum(w.sum(1,keepdims=True),1e-9)
    return np.log(np.nansum(w*np.nan_to_num(Rr),1))
def sist(e,cel,mask):
    ids,inv=np.unique(cel[mask],return_inverse=True); n=np.bincount(inv); me=np.bincount(inv,e[mask])/n
    return np.sqrt(np.sum(n*me**2)/n.sum())
cel=j*1000+i
tr=ok&~cass; te=ok&cass
def valuta(X, nome):
    X1=np.column_stack([np.ones(tr.sum())]+[x[tr] for x in X]); c=np.linalg.lstsq(X1,y[tr],rcond=None)[0]
    p=c[0]+sum(ci*x[te] for ci,x in zip(c[1:],X))
    e=(mod[te]*np.exp(p)-vero[te])/0.5144
    v=dist[te]<3
    print(f"{nome:34s} coeff {np.round(c,2).tolist()} rmse {np.sqrt(np.mean(e**2)):.3f} sist {sist(e,cel[te],np.ones(te.sum(),bool)):.3f} | <3km rmse {np.sqrt(np.mean(e[v]**2)):.3f} sist {sist(e,cel[te],v):.3f}")
valuta([], "modello")
valuta([lL], "lens")
for sd in [0,15,30]:
    for ss in [0,1,2]:
        lc=lc_var(sd,ss)
        valuta([lc], f"cfd sig_dir {sd} sig_spazio {ss}")
        valuta([lL,lc], f"lens+cfd sig_dir {sd} sig_spazio {ss}")
