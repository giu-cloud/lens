"""Il livello del radar sul bordo del riquadro: quanto è lontano da 1 dove LensWind raccorda al modello (1,5 miglia)."""
import sys, subprocess, numpy as np, os, datetime as dt
sys.path.insert(0,"/home/claude/lens/banco/cfd")
from riferimenti import AREE, leggi_rif, CP, MAPPE, L, D
km=sys.argv[1]
for nome in sys.argv[2:]:
    c=AREE[nome]; T=np.load(f"{L}/{c['tab']}.npz",mmap_mode="r"); sel=np.asarray(T["area"])==c["idx"]
    t=np.asarray(T["t"][sel]); tt=np.unique(t[(t//(7*86400))%5==2])
    tt=tt[np.random.default_rng(3).permutation(len(tt))[:24]]
    sc=[]
    for x in tt:
        d=dt.datetime.fromtimestamp(int(x),dt.timezone.utc); run=d.replace(hour=0 if d.hour<12 else 12,minute=0,second=0)
        pw=f"{D}/ecmwf_ampio/ecmwf_{run:%Y%m%d%H}.grb2"; pp=f"{D}/stabw/stabw_{run:%Y%m%d%H}.grb2"
        if os.path.exists(pw): sc+= [pw,pp,str(int(x)*1000),"nan","nan"]
    usc=f"/tmp/rif/bordo_{nome}.bin"
    r=subprocess.run(["java","-Xmx3g","-Dshift.lens.relief=0.25","-Dlens.soloFast=true","-cp",CP,"RiferimentiKt",c["atl"],c["lib"],MAPPE,km if nome!="messina" else "5",usc]+sc,capture_output=True,text=True)
    box,nx,ny,out=leggi_rif(usc); os.remove(usc)
    S,N,O,E=box; dy=(N-S)/ny*60; dx=(E-O)/nx*60*np.cos(np.radians((S+N)/2))
    jj,ii=np.meshgrid(np.arange(ny),np.arange(nx),indexing="ij")
    bordo=np.minimum(np.minimum(jj,ny-1-jj)*dy,np.minimum(ii,nx-1-ii)*dx)   # miglia dal bordo
    fb=[];fi=[];salto=[]
    for t2,fr,n,F in out:
        mare=(F[9]>0.5)&np.isfinite(F[1])&(F[1]>0.3)
        liv=F[6]; tot=np.hypot(F[2],F[3])*F[6]/F[1]
        a=mare&(bordo<0.3); b=mare&(bordo>3)
        if a.sum()<20: continue
        fb.append(np.abs(tot[a]-1)); fi.append(np.abs(tot[b]-1)); salto.append(np.percentile(np.abs(tot[a]-1),[50,90,99]))
    fb=np.concatenate(fb)
    print(f"{nome} km {km}: sul bordo |C/modello - 1|: mediana {np.median(fb):.3f}  p90 {np.percentile(fb,90):.3f}  p99 {np.percentile(fb,99):.3f}; quota oltre 15%: {np.mean(fb>0.15)*100:.1f}%  oltre 25%: {np.mean(fb>0.25)*100:.1f}%")
