"""Figura del Pellegrino: CFD, atlante di Lens con il rilievo acceso e la carta di Pietro, in rapporto al vento al largo."""
import os
import numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'banco/cfd'); sys.path.insert(0,'banco/strumenti')
from controlla_atlante import leggi
from taratura_radar import mare_su
A=leggi('/home/claude/dati/shift/atlante_palermo.lna'); s,n,w,e=A['box']; ny,nx=A['ny'],A['nx']
lat=n-(np.arange(ny)+0.5)*(n-s)/ny; lon=w+(np.arange(nx)+0.5)*(e-w)/nx
mare=mare_su('/home/claude/dati/shift/atlante_palermo.lna',lat,lon)
zC=np.load('/home/claude/dati/cfd/palermo/cfd_palermo_281.2.npz')
d=278; u0,v0=-np.sin(np.radians(d)),-np.cos(np.radians(d))
fr=[float(x) for x in A['fr']]; k=fr.index(min(fr,key=lambda x:abs(x-0.7)))
ue,ve,un,vn=A['campi'][k]; Ua=u0*ue+v0*un; Va=u0*ve+v0*vn
b=np.load('/tmp/pietro_bin.npy'); V=b[0,16].reshape(ny,nx)/2.0; Dd=b[2,16].reshape(ny,nx)*360/128
off=np.median(V[(lat[:,None]>38.23)&(lon[None,:]<13.36)&mare])
Up=-V/off*np.sin(np.radians(Dd)); Vp=-V/off*np.cos(np.radians(Dd))
fig,ax=plt.subplots(1,3,figsize=(18,6.5))
for a,(tit,U,Vv) in zip(ax,[("CFD (vento da 281°)",zC['u'],zC['v']),("Atlante di Lens con il rilievo acceso (278°, Fr 0,7)",Ua,Va),("Carta di Pietro, gio 26 ore 03 (diviso il vento al largo)",Up,Vp)]):
    S=np.where(mare,np.hypot(U,Vv),np.nan)
    im=a.pcolormesh(lon,lat,S,vmin=0.3,vmax=1.3,cmap='viridis',shading='auto')
    a.contour(lon,lat,mare.astype(float),[0.5],colors='k',linewidths=0.7)
    kk=4; sel=np.zeros_like(mare); sel[::kk,::kk]=True; sel&=mare
    LA,LO=np.meshgrid(lat,lon,indexing='ij')
    a.quiver(LO[sel],LA[sel],(U/np.hypot(U,Vv))[sel],(Vv/np.hypot(U,Vv))[sel],scale=45,width=0.0025)
    a.set_xlim(13.33,13.43); a.set_ylim(38.13,38.23); a.set_title(tit,fontsize=10); a.set_aspect(1/np.cos(np.radians(38.18)))
    for nome,(la_,lo_) in {"P. del Rotolo":(38.192,13.362),"Vergine Maria":(38.170,13.368),"Arenella":(38.150,13.373)}.items():
        a.plot(lo_,la_,'r.',ms=5); a.text(lo_+0.002,la_,nome,fontsize=7,color='r')
plt.colorbar(im,ax=ax,shrink=0.7,label='vento / vento al largo')
plt.savefig('/home/claude/lens/banco/figure/pellegrino_cfd_atlante_carta.png',dpi=90,bbox_inches='tight')
