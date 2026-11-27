import numpy as np
L="/home/claude/lavoro/cfd_rif"
AREE=["messina","palermo","eolie","cagliari","monaco","corsica-nord","dalmazia"]
KS=["0","0.25","0.5","0.75",""]
KV=np.array([0,.25,.5,.75,1.])
def carica(a):
    R=[np.load(f"{L}/rif_{a}"+(f"_r{k}" if k else "")+".npz") for k in KS]
    n=[len(r["t"]) for r in R]
    assert len(set(n))==1,(a,n)
    for r in R[1:]: assert np.array_equal(r["t"],R[0]["t"]) and np.array_equal(r["j"],R[0]["j"])
    d={k:np.asarray(R[-1][k]) for k in ("scena","t","j","i","lat","lon","dist","vero","mod","Fr","dirM")}
    d["r"]=np.stack([np.asarray(r["rFast"]) for r in R],1)
    d["ok"]=np.isfinite(d["vero"])&(d["vero"]>0.5)&(d["mod"]>2.0)&np.isfinite(d["r"]).all(1)
    d["cel"]=d["j"]*1000+d["i"]
    return d
def sist(cel,e,mk):
    ids,inv=np.unique(cel[mk],return_inverse=True); n=np.bincount(inv); me=np.bincount(inv,e[mk])/n
    return float(np.sqrt(np.sum(n*me**2)/n.sum()))
