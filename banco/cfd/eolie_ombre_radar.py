import numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'banco/cfd'); sys.path.insert(0,'banco/strumenti')
from riva_radar import risposta_su_radar, AREE14
from controlla_atlante import leggi
T=np.load('/home/claude/lavoro/tabella6.npz', mmap_mode='r'); V=np.load('/home/claude/lavoro/radar_vero_tabella6.npy', mmap_mode='r')
m=np.asarray(T['area'])==AREE14.index('eolie')
C=np.load('/home/claude/lavoro/base3/campi_eolie.npz'); rlat=C['lat'].astype(float); rlon=C['lon'].astype(float)
R=risposta_su_radar('/home/claude/dati/shift/atlante_eolie.lna',rlat,rlon)
j=T['j'][m].astype(int); i=T['i'][m].astype(int); vero=V[m]; mod=T['mod'][m]; q=T['q_lens054'][m]; S0=T['S0'][m]; D0=T['D0'][m]; dT=T['dT925'][m]; sc=T['scena'][m]; dist=T['dist'][m]
x=(D0%360)/22.5; k0=np.floor(x).astype(int)%16; w=x-np.floor(x); k1=(k0+1)%16
ra=(1-w)*R[k0,j,i]+w*R[k1,j,i]
A=leggi('/home/claude/dati/shift/atlante_eolie.lna'); s,n,ww,e=A['box']; mm=A['m']; mare=A['mare'][mm:A['NY']-mm,mm:A['NX']-mm]
la=n-(np.arange(A['ny'])+0.5)*(n-s)/A['ny']; lo=ww+(np.arange(A['nx'])+0.5)*(e-ww)/A['nx']
setts=[("ponente-maestrale 270-330",(270,330)),("tramontana-grecale 330-40",(330,40)),("scirocco-libeccio 120-240",(120,240))]
fig,ax=plt.subplots(len(setts),4,figsize=(20,5*len(setts)))
ok=np.isfinite(vero)&(mod>2)&(q>0)&np.isfinite(ra)&(S0>7)
for r,(nome,(a,b)) in enumerate(setts):
    sel=ok&(((D0>=a)&(D0<b)) if a<b else ((D0>=a)|(D0<b)))
    stab=sel&(dT>2.0)
    for c,(tit,x_,f) in enumerate([("radar / largo",vero/S0,sel),("Lens 0.54 / largo",mod*q/S0,sel),("rilievo pieno / largo",mod*ra/S0,sel),("radar / largo, aria stabile",vero/S0,stab)]):
        G=np.full((len(rlat),len(rlon)),np.nan); N=np.zeros_like(G); Ssum=np.zeros_like(G)
        np.add.at(N,(j[f],i[f]),1); np.add.at(Ssum,(j[f],i[f]),x_[f])
        G=np.where(N>=3,Ssum/np.maximum(N,1),np.nan)
        im=ax[r,c].pcolormesh(rlon,rlat,G,cmap='RdYlBu',vmin=0.5,vmax=1.3,shading='auto')
        ax[r,c].contour(lo,la,mare.astype(float),[0.5],colors='k',linewidths=0.5)
        ax[r,c].set_title(f"{nome}: {tit} ({len(np.unique(sc[f]))} pass.)",fontsize=9)
        ax[r,c].set_aspect(1/np.cos(np.radians(38.5)))
plt.colorbar(im,ax=ax.ravel().tolist(),shrink=0.5,label='rapporto con il vento al largo')
plt.savefig('/home/claude/lens/banco/figure/eolie_ombre_radar.png',dpi=65,bbox_inches='tight')
# numeri: sottovento entro 1,5-4 km, per stabilità
terra=T['terra'][m]; ang=T['ang'][m]
for lab,f in (("tutte",ok),("stabile dT925>2",ok&(dT>2)),("rimescolata dT925<0",ok&(dT<0))):
    for lo_,hi_ in ((1.5,2.5),(2.5,4),(4,8)):
        g=f&(terra<3)&(ang>4)&(dist>=lo_)&(dist<hi_)
        if g.sum()<50: continue
        print(f"{lab:22s} {lo_}-{hi_} km sottovento n {g.sum():5d}: radar/largo {np.median(vero[g]/S0[g]):.2f} | Lens0.54/largo {np.median(mod[g]*q[g]/S0[g]):.2f} | modello/largo {np.median(mod[g]/S0[g]):.2f} | rilievo pieno/largo {np.median(mod[g]*ra[g]/S0[g]):.2f}")
