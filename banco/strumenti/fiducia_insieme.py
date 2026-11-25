"""La fiducia dell'insieme: quale misura (dispersione di ECMWF, di ICON-2I, media 0,6/0,4, disaccordo dei modelli) segue meglio l'errore dell'insieme corretto alle stazioni."""
import json,glob,numpy as np
R=[]
for f in glob.glob('/home/claude/lavoro/stazioni_prova_r8b/*.json'): R+=json.load(open(f))['righe']
g=lambda k: np.array([r.get(k,np.nan) for r in R],float)
o,ec,ic,de,di,S0=g('staz_mare'),g('ecmwf_c'),g('icon_c'),g('disp'),g('disp_i'),g('S0')
ens=0.6*ec+0.4*ic
ok=np.isfinite(ens)&np.isfinite(de)&np.isfinite(di)&(S0>=2)&(o>1)
o,ec,ic,de,di,ens=o[ok],ec[ok],ic[ok],de[ok],di[ok],ens[ok]
rel=np.abs(ens-o)/np.maximum(o,2.0)
disacc=np.abs(ec-ic)/np.maximum(ens,2.0)
cand={"dispersione di ECMWF":de,"dispersione di ICON-2I":di,"media 0,6/0,4":0.6*de+0.4*di,
      "media + disaccordo dei modelli":0.6*de+0.4*di+disacc,"solo disaccordo dei modelli":disacc}
from scipy.stats import spearmanr
print("confronti",ok.sum())
for k,x in cand.items():
    q=np.quantile(x,[0.2,0.4,0.6,0.8]); c=np.digitize(x,q)
    med=[np.median(rel[c==i]) for i in range(5)]
    print(f"{k:34s} Spearman {spearmanr(x,rel).correlation:.3f}  errore relativo per quinto: "+" ".join(f"{m:.2f}" for m in med)+f"  rapporto ultimo/primo {med[4]/med[0]:.2f}")
# soglie: con la media 0,6/0,4 + disaccordo, quali soglie danno le stesse tre classi della carta?
x=0.6*de+0.4*di+disacc
print("quantili 20/40/60/80 della media+disaccordo:", np.quantile(x,[0.2,0.4,0.6,0.8]).round(3), " della media 0,6/0,4:", np.quantile(0.6*de+0.4*di,[0.2,0.4,0.6,0.8]).round(3))
for lo,hi in ((0,0.22),(0.22,0.33),(0.33,9)):
    m=(0.6*de+0.4*di>=lo)&(0.6*de+0.4*di<hi); print(f"  media 0,6/0,4 in [{lo},{hi}): n={m.sum()} errore relativo mediano {np.median(rel[m]):.2f}")
