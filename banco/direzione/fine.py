import sys, numpy as np
from analisi import dire, diff
L="/home/claude/lavoro/dir"
print("stazione: n | errore ass. medio sul mare vicino: modello, Lens atlante grande, Lens atlante a 100 m | nella cella della stazione: grande, 100 m")
for a in sys.argv[1:]:
    G=np.load(f"{L}/dir_{a}.npz"); F=np.load(f"{L}/dir_{a}_fine.npz")
    kg={(s,x[0]):x for s,x in zip(G["st"],G["X"])}
    for st in np.unique(F["st"]):
        rows=[(kg[(s,x[0])],x) for s,x in zip(F["st"],F["X"]) if s==st and (s,x[0]) in kg]
        if len(rows)<60: continue
        g=np.array([r[0] for r in rows]); f=np.array([r[1] for r in rows])
        ok=(np.hypot(g[:,15],g[:,16])>=5)&(g[:,1]>=3)&np.isfinite(f[:,19])&np.isfinite(g[:,19])
        ds=g[ok,2]
        e=lambda X,o: np.abs(diff(ds,dire(X[ok,o],X[ok,o+1]))).mean()
        print(f"{a:13s} {st:16s} {ok.sum():4d} | {e(g,15):5.1f}  {e(g,19):5.1f}  {e(f,19):5.1f} | {e(g,9):5.1f}  {e(f,9):5.1f}")
