"""Riassunto di raffiche_prova.py: la raffica delle stazioni segue il fattore di Lens 2? Pendenze dentro stazione e ottante, ed errore libero delle varianti."""
import json,glob,numpy as np
KN=1.943844
R=[]
for f in glob.glob('/home/claude/lavoro/raffiche/*.json'): R+=json.load(open(f))
A={k:np.array([r[k] for r in R]) for k in R[0]}
ok=(A['mod']>=3)&(A['S0']>=2)&(A['staz']>0.5)&(A['staz_raff']>0.5)&(A['mod_raff']>1)
A={k:v[ok] for k,v in A.items()}
grp=np.array([f"{s}_{o}" for s,o in zip(A['st'],A['ott'])])
print("confronti", ok.sum(), "gruppi stazione×ottante", len(set(grp)))
def demean(x):
    y=x.copy()
    for g in set(grp):
        m=grp==g; y[m]-=x[m].mean()
    return y
def pend(y,x):
    yd,xd=demean(y),demean(x)
    # pendenze separate sotto e sopra 1
    Xm=np.minimum(x,0); Xp=np.maximum(x,0)
    X=np.stack([demean(Xm),demean(Xp)],1)
    b,res,_,_=np.linalg.lstsq(X,yd,rcond=None)
    e=yd-X@b; s2=e@e/(len(yd)-2); cov=s2*np.linalg.inv(X.T@X); se=np.sqrt(np.diag(cov))
    b1=np.sum(xd*yd)/np.sum(xd*xd); se1=np.sqrt(np.sum((yd-b1*xd)**2)/(len(yd)-1)/np.sum(xd*xd))
    return b1,se1,b,se
lg=np.log(A['G']); lgr=np.log(A['G_raff'])
for nome,y,x in (("forza media",np.log(A['staz']/A['mod']),lg),("raffica",np.log(A['staz_raff']/A['mod_raff']),lgr)):
    b1,se1,b,se=pend(y,x)
    print(f"{nome:12s} pendenza su log(fattore): {b1:.2f} ± {se1:.2f}   sotto 1: {b[0]:.2f} ± {se[0]:.2f}   sopra 1: {b[1]:.2f} ± {se[1]:.2f}")
print("fattore: quota sotto 0,95", np.mean(A['G']<0.95).round(2), "sopra 1,05", np.mean(A['G']>1.05).round(2), "min/max", A['G'].min().round(2), A['G'].max().round(2))
def libero(o,p):
    num=0;den=0
    for g in set(grp):
        m=grp==g
        if m.sum()<10: continue
        k=np.sum(o[m]*p[m])/np.sum(p[m]*p[m]); num+=np.sum((k*p[m]-o[m])**2); den+=m.sum()
    return np.sqrt(num/den)*KN
o=A['staz_raff']; g0=A['mod_raff']; G=A['G_raff']
var={"raffica del modello":g0,"× fattore":g0*G,"× fattore solo sopra 1":g0*np.maximum(G,1),"× radice del fattore":g0*np.sqrt(G),
     "× fattore sopra 1, metà sotto":g0*np.where(G>=1,G,(1+G)/2)}
print("\nraffica, errore quadratico libero (nodi), tutte / dove il fattore è sotto 0,95 / sopra 1,05:")
for k,p in var.items():
    m1=G<0.95; m2=G>1.05
    def lib_m(m):
        num=0;den=0
        for g in set(grp[m]):
            mm=m&(grp==g)
            if mm.sum()<5: continue
            # scala presa su tutta la serie del gruppo, poi errore solo sul sottoinsieme
            mg=grp==g; kk=np.sum(o[mg]*p[mg])/np.sum(p[mg]*p[mg]); num+=np.sum((kk*p[mm]-o[mm])**2); den+=mm.sum()
        return np.sqrt(num/den)*KN if den else float('nan')
    print(f"  {k:32s} {libero(o,p):.3f}   {lib_m(m1):.3f} (n={m1.sum()})   {lib_m(m2):.3f} (n={m2.sum()})")
print("\nforza media, errore quadratico libero: modello", round(libero(A['staz'],A['mod']),3), "corretto", round(libero(A['staz'],A['mod']*A['G']),3))
for s in sorted(set(A['nome'])):
    m=A['nome']==s
    print(f"  {s:16s} n={m.sum():4d} fattore medio {A['G'][m].mean():.2f}  raffica staz/mod {np.median(o[m]/g0[m]):.2f}  forza staz/mod {np.median(A['staz'][m]/A['mod'][m]):.2f}")
