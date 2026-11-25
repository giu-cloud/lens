"""Carte medie per settore di vento: radar/modello, CFD/largo, Lens/modello, radar/largo (uscita di taratura_radar.py --npz).

  SETTORI="[(225,260),(260,300)]" python3 mappe_settori.py <righe .npz> <uscita .png> <terreno .npz o atlante .lna> <un cfd_<sito>_XXX.X.npz>
"""
import numpy as np, matplotlib, sys; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'banco/cfd'); from taratura_radar import mare_su
f=sys.argv[1]; out=sys.argv[2]; mfonte=sys.argv[3]; cfdz=sys.argv[4]
z=np.load(f)
ok=z['ok']; lL=z['lL']; vero=z['vero']; mod=z['mod']; D0=z['D0']; S0=z['S0']; sc=z['scena']; j=z['j']; i=z['i']; lC=z['lC']
rlat=z['rlat']; rlon=z['rlon']
c0=np.load(cfdz); la, lo=c0['lat'], c0['lon']; mare=mare_su(mfonte, la, lo)
import os
setts=eval(os.environ.get('SETTORI','[(225,260),(260,300),(300,340),(0,60)]'))
fig,ax=plt.subplots(len(setts),4,figsize=(18,4.2*len(setts)))
for r,(a,b) in enumerate(setts):
    k=ok&(D0>=a)&(D0<b)&(S0>4)
    for c,(nome,x) in enumerate([("radar/modello",np.log(vero/mod)),("CFD/largo",lC),("Lens/modello",lL),("radar/largo",np.log(vero/S0))]):
        G=np.full((len(rlat),len(rlon)),np.nan); N=np.zeros_like(G)
        np.add.at(N,(j[k],i[k]),1); S=np.zeros_like(G); np.add.at(S,(j[k],i[k]),x[k])
        G=np.where(N>=4,S/np.maximum(N,1),np.nan); G=G-np.nanmedian(G)
        im=ax[r,c].pcolormesh(rlon,rlat,G,cmap='RdBu_r',vmin=-0.3,vmax=0.3)
        ax[r,c].contour(lo,la,mare.astype(float),[0.5],colors='k',linewidths=0.5)
        ax[r,c].set_title(f"{a}-{b}°: {nome} ({len(np.unique(sc[k]))} pass.)",fontsize=9)
        ax[r,c].set_xlim(lo.min(),lo.max()); ax[r,c].set_ylim(la.min(),la.max())
plt.colorbar(im,ax=ax.ravel().tolist(),shrink=0.5,label='log, scarto dalla mediana')
plt.savefig(out,dpi=70)
