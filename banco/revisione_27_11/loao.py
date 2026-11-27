from carica import *
import itertools
D={a:carica(a) for a in AREE}
def pred(d,k):  # k per riga, interpola tra le manopole salvate
    k=np.clip(k,0,1); x=k/0.25; i0=np.minimum(x.astype(int),3); f=x-i0
    r=d["r"][np.arange(len(k)),i0]*(1-f)+d["r"][np.arange(len(k)),i0+1]*f
    return d["mod"]*r
def regola(fr,p):
    kl,kh,f1,f2=p
    return kl+(kh-kl)*np.clip((fr-f1)/max(f2-f1,1e-6),0,1)
def punteggi(d,k):
    e=(pred(d,k)-d["vero"])/0.5144; ok=d["ok"]; m3=ok&(d["dist"]<3)
    return np.sqrt(np.mean(e[m3]**2)), sist(d["cel"],e,m3), np.sqrt(np.mean(e[ok]**2)), sist(d["cel"],e,ok)
def costo(aree,p):
    s=0
    for a in aree:
        d=D[a]; q=punteggi(d,regola(d["Fr"],p)); b=BASE[a]
        s+= (q[0]/b[0])**2+(q[1]/b[1])**2+(q[2]/b[2])**2
    return s
BASE={a:punteggi(D[a],np.full(len(D[a]["t"]),0.25)) for a in AREE}
griglia=[(kl,kh,f1,f2) for kl in (0.25,0.5,0.75,1.0) for kh in (0.0,0.1,0.25) for f1 in (0.15,0.2,0.3,0.4) for f2 in (0.3,0.4,0.5,0.7,1.0) if f2>f1]
print("area          0.25: rms3 sist3 rms sist | regola tarata sulle altre sei")
tot=[]
for a in AREE:
    altre=[x for x in AREE if x!=a]
    best=min(griglia,key=lambda p:costo(altre,p))
    q=punteggi(D[a],regola(D[a]["Fr"],best)); b=BASE[a]; z=punteggi(D[a],np.zeros(len(D[a]["t"])))
    print(f"{a:13s} {b[0]:.2f} {b[1]:.2f} {b[2]:.2f} {b[3]:.2f} | {q[0]:.2f} {q[1]:.2f} {q[2]:.2f} {q[3]:.2f}  p={best} | k=0: {z[0]:.2f} {z[1]:.2f} {z[2]:.2f} {z[3]:.2f}")
best=min(griglia,key=lambda p:costo(AREE,p)); print("tutte:",best, costo(AREE,best), "0.25:",costo(AREE,(0.25,0.25,0.2,0.3)), "0:",costo(AREE,(0,0,0.2,0.3)))
for p in sorted(griglia,key=lambda p:costo(AREE,p))[:8]: print(p, round(costo(AREE,p),3))
