import sys, numpy as np
from carica import AREE, sist
L="/home/claude/lavoro/cfd_rif"
def leggi(a,suf):
    R=np.load(f"{L}/rif_{a}{suf}.npz"); return {k:np.asarray(R[k]) for k in R}
def tabella(sufs, base="_r0.25", boot=True):
    rng=np.random.default_rng(1)
    print("sistematico entro 3 km / quadratico entro 3 km / quadratico tutto ; differenza dalla base con intervallo 90% (sist3, rms tutto)")
    for a in AREE:
        B=leggi(a,base); out=[]
        for s in [base]+sufs:
            try: R=leggi(a,s)
            except FileNotFoundError: out.append(f"{s}: -"); continue
            assert np.array_equal(R["t"],B["t"])
            ok=np.isfinite(B["vero"])&(B["vero"]>0.5)&(B["mod"]>2)&np.isfinite(R["rFast"])&np.isfinite(B["rFast"])
            m3=ok&(B["dist"]<3); cel=B["j"]*1000+B["i"]
            e=(R["mod"]*R["rFast"]-R["vero"])/0.5144; eb=(B["mod"]*B["rFast"]-B["vero"])/0.5144
            txt=f"{s}: {sist(cel,e,m3):.2f}/{np.sqrt(np.mean(e[m3]**2)):.2f}/{np.sqrt(np.mean(e[ok]**2)):.2f}"
            if s!=base and boot:
                sc=np.unique(B["scena"][ok]); idx={x:np.where(ok&(B["scena"]==x))[0] for x in sc}
                d1=[];d2=[]
                for _ in range(200):
                    pick=np.concatenate([idx[x] for x in rng.choice(sc,len(sc))])
                    mk=np.zeros(len(ok),bool); 
                    p3=pick[B["dist"][pick]<3]
                    def s_(ee,pp):
                        ids,inv=np.unique(cel[pp],return_inverse=True); n=np.bincount(inv); me=np.bincount(inv,ee[pp])/n
                        return np.sqrt(np.sum(n*me**2)/n.sum())
                    d1.append(s_(e,p3)-s_(eb,p3)); d2.append(np.sqrt(np.mean(e[pick]**2))-np.sqrt(np.mean(eb[pick]**2)))
                q1=np.percentile(d1,[5,95]); q2=np.percentile(d2,[5,95])
                txt+=f" [{q1[0]:+.2f}..{q1[1]:+.2f} | {q2[0]:+.2f}..{q2[1]:+.2f}]"
            out.append(txt)
        print(f"{a:13s} "+"  ".join(out))
if __name__=="__main__":
    tabella(sys.argv[2:], sys.argv[1])
