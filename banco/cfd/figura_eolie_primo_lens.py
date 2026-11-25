import numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/tmp/verifica'); sys.path.insert(0,'/home/claude/lens/banco/strumenti')
from leggi_campi import leggi
from scipy.ndimage import distance_transform_edt
from controlla_atlante import leggi as leggi_lna
venti=[(315,15),(270,15)]
la1,lo1,F1=leggi('/tmp/verifica/eolie340.bin',2)
la2,lo2,F2=leggi('/tmp/verifica/eolie100.bin',2)
def mare_di(p):
    A=leggi_lna(p); m=A['m']; return A['mare'][m:A['NY']-m,m:A['NX']-m], A['dx']
M1,dx1=mare_di('/home/claude/dati/shift/atlante_eolie.lna'); M2,dx2=mare_di('/tmp/verifica/atlante_eolie-centro-100.lna')
d2=distance_transform_edt(M2)*dx2/1000
def peso(d,d1=1.0,dd2=2.5):
    x=np.clip((d-d1)/(dd2-d1),0,1); return 0.5+0.5*np.cos(np.pi*x)
zoom=(38.36,38.60,14.80,15.04)
fig,ax=plt.subplots(len(venti),4,figsize=(24,6.4*len(venti)))
for r,(da,kn) in enumerate(venti):
    (u0,v0),(u1,v1)=F2[r]; (a0,b0),(a1,b1)=F1[r]
    s_oggi=np.hypot(u0,v0)*1.944; s_primo=np.hypot(u1,v1)*1.944
    s_riva=s_oggi*(np.maximum(s_primo,0.05)/np.maximum(s_oggi,0.05))**peso(d2)
    pannelli=[("Lens di oggi (0.56), 100 m",la2,lo2,s_oggi,(u0,v0),M2),
              ("Primo Lens com'era: 340 m",la1,lo1,np.hypot(a1,b1)*1.944,(a1,b1),M1),
              ("Primo Lens su griglia di 100 m",la2,lo2,s_primo,(u1,v1),M2),
              ("Primo Lens solo vicino alla riva (pieno entro 1 km, sfuma a 2,5)",la2,lo2,s_riva,(u0,v0),M2)]
    for c,(tit,la,lo,S,(U,V),M) in enumerate(pannelli):
        a=ax[r,c]
        im=a.pcolormesh(lo,la,np.where(M,S,np.nan),cmap='turbo',vmin=0,vmax=24,shading='auto')
        a.contour(lo,la,M.astype(float),[0.5],colors='k',linewidths=0.6)
        k=max(1,int(round(1.2/((la[0]-la[1])*111.19))))
        LA,LO=np.meshgrid(la,lo,indexing='ij'); sel=np.zeros_like(M); sel[::k,::k]=True; sel&=M
        sp=np.hypot(U,V)+1e-9
        a.quiver(LO[sel],LA[sel],(U/sp)[sel],(V/sp)[sel],scale=45,width=0.002)
        a.set_xlim(zoom[2],zoom[3]); a.set_ylim(zoom[0],zoom[1]); a.set_aspect(1/np.cos(np.radians(38.5)))
        a.set_title(f"da {da}°, {kn} nodi, Fr 0,7: {tit}",fontsize=9)
plt.colorbar(im,ax=ax.ravel().tolist(),shrink=0.6,label='nodi')
plt.savefig('/home/claude/lens/banco/figure/eolie_primo_lens_fine.png',dpi=70,bbox_inches='tight')
# numeri in punti noti
punti={"Bocche di Vulcano (canale)":(38.431,14.952),"sotto Lipari, costa est (Canneto)":(38.495,14.965),"sotto Lipari, costa sud-est (Marina Corta, 500 m)":(38.462,14.968),"largo a ovest":(38.50,14.82)}
for r,(da,kn) in enumerate(venti):
    (u0,v0),(u1,v1)=F2[r]
    for nome,(p,q) in punti.items():
        j=np.argmin(abs(la2-p)); i=np.argmin(abs(lo2-q))
        print(f"{da}° {nome:48s} mare {M2[j,i]} oggi {np.hypot(u0[j,i],v0[j,i])*1.944:5.1f} kn | primo Lens 100 m {np.hypot(u1[j,i],v1[j,i])*1.944:5.1f} kn | dist {d2[j,i]:.2f} km")
