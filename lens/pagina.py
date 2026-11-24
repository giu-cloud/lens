"""Pagina della carta del vento sui golfi: un solo file HTML, con le ore da scorrere.
Lettura pensata per il telefono: colori a fasce di forza Beaufort (in nodi), banderuole su una griglia fissa,
tocco sulla carta per leggere intensita' e direzione nel punto, vento mostrato sul mare (a terra a richiesta)."""
import json

CSS = """
:root{color-scheme:light dark;--bg:#f4f3ef;--card:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--line:#e3e1da}
@media (prefers-color-scheme:dark){:root{--bg:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--line:#2e2e2b}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:820px;margin:0 auto;padding:14px 12px 40px}h1{font-size:1.35rem;margin:.2em 0 0}
.sub{color:var(--ink2);font-size:.88rem;margin:0 0 10px}
.card{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:12px;margin:0 0 12px}
.riga{display:flex;gap:6px;overflow-x:auto;padding-bottom:8px}
button{font:inherit;font-size:.85rem;border:1px solid var(--line);background:var(--card);color:var(--ink);border-radius:999px;padding:6px 11px;cursor:pointer;white-space:nowrap}
button[aria-pressed=true]{background:var(--ink);color:var(--card);border-color:var(--ink)}
canvas{display:block;width:100%;border-radius:10px;touch-action:manipulation;background:#c6def0}
.tempo{display:flex;align-items:center;gap:8px;margin:10px 0 2px}.tempo input{flex:1;min-width:0;accent-color:#2a78d6}
.tl{font-weight:650;font-size:1.05rem;font-variant-numeric:tabular-nums}.reg{color:var(--ink2);font-size:.86rem;margin:2px 0 8px}
.bande{display:grid;grid-template-columns:repeat(8,1fr);gap:2px;margin:8px 0 2px;text-align:center;font-size:.72rem;line-height:1.2;font-variant-numeric:tabular-nums}
.bande i{display:block;height:12px;border-radius:3px;margin-bottom:2px}.bande small{display:block;color:var(--ink2);font-size:.66rem}
.scala{height:10px;border-radius:5px;margin:8px 0 2px}.tacche{display:flex;justify-content:space-between;font-size:.72rem;color:var(--ink2)}
.info{font-size:.92rem;margin:8px 0 0;min-height:2.8em}.note{color:var(--ink2);font-size:.85rem}.note li{margin:.35em 0}
.note svg{vertical-align:middle}
"""

JS = r"""
(async function(){
const D=__DATI__,Q=id=>document.getElementById(id);
function dec(b){const s=atob(b),a=new Uint8Array(s.length);for(let i=0;i<s.length;i++)a[i]=s.charCodeAt(i);return a}
let P;
try{P=JSON.parse(await new Response(new Blob([dec(D.gz)]).stream().pipeThrough(new DecompressionStream('gzip'))).text())}
catch(e){Q('info').textContent='Questo browser è troppo vecchio per aprire i dati della carta: aggiornalo e riprova.';return}
const ORE=P.ore,nx=D.nx,ny=D.ny,KB=D.kb,N=nx*ny;
const quota=new Uint16Array(dec(P.quota).buffer);
const MAR=new Uint8Array(N);{const m=dec(P.mare);for(let k=0;k<N;k++)MAR[k]=(m[k>>3]>>(7-(k&7)))&1}
function bl(a,i,j){const gx=(i+D.m+.5)/D.NX*(D.gx-1),gy=(j+D.m+.5)/D.NY*(D.gy-1),x0=Math.min(D.gx-2,Math.floor(gx)),y0=Math.min(D.gy-2,Math.floor(gy)),wx=gx-x0,wy=gy-y0,q=y0*D.gx+x0;
 return(a[q]*(1-wx)+a[q+1]*wx)*(1-wy)+(a[q+D.gx]*(1-wx)+a[q+D.gx+1]*wx)*wy}
function ora(n){const o=ORE[n];if(o.S)return o;
 const V=dec(o.v),Dd=dec(o.d),Rb=dec(o.r),G=o.g?dec(o.g):V;
 o.S=new Float32Array(N);o.G=new Float32Array(N);o.U=new Float32Array(N);o.W=new Float32Array(N);o.R=new Uint8Array(N);
 o.MS=new Float32Array(N);o.MU=new Float32Array(N);o.MW=new Float32Array(N);
 const mu=o.mv.map((s,q)=>-s*Math.sin(o.md[q]*Math.PI/180)),mw=o.mv.map((s,q)=>-s*Math.cos(o.md[q]*Math.PI/180));
 for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const k=j*nx+i,s=V[k]/4,a=Dd[k]*Math.PI/128;
  o.S[k]=s;o.G[k]=G[k]/4;o.U[k]=-s*Math.sin(a);o.W[k]=-s*Math.cos(a);o.R[k]=(Rb[k>>3]>>(7-(k&7)))&1;
  o.MS[k]=bl(o.mv,i,j);o.MU[k]=bl(mu,i,j);o.MW[k]=bl(mw,i,j)}
 return o}
// scala del vento di Shift
const WS=[[0,[127,178,229]],[5,[76,175,80]],[10,[192,202,51]],[15,[255,179,0]],[20,[255,112,67]],[25,[229,57,53]],[32,[142,36,170]]];
function wcol(v){if(v<=0)return WS[0][1];for(let i=1;i<WS.length;i++)if(v<=WS[i][0]){const a=WS[i-1],b=WS[i],t=(v-a[0])/(b[0]-a[0]);return[a[1][0]+(b[1][0]-a[1][0])*t,a[1][1]+(b[1][1]-a[1][1])*t,a[1][2]+(b[1][2]-a[1][2])*t]}return WS[WS.length-1][1]}
const rgb=c=>'rgb('+Math.round(c[0])+','+Math.round(c[1])+','+Math.round(c[2])+')';
const LIM=[3.5,6.5,10.5,16.5,21.5,27.5,33.5,40.5],BF=['F0–1','F2','F3','F4','F5','F6','F7','F8','F9+'];
function forza(v){let b=0;while(b<LIM.length&&v>=LIM[b])b++;return BF[b]}
const num=(x,n)=>x.toFixed(n||0).replace('.',',').replace('-','−');
const ROSA=['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSO','SO','OSO','O','ONO','NO','NNO'];
const NOMI=['tramontana','grecale','levante','scirocco','ostro','libeccio','ponente','maestrale'];
const rosa=d=>ROSA[Math.round(d/22.5)%16],nome=d=>NOMI[Math.round(d/45)%8];
let h=0,strato='ril',campitura=false,barbe=false,vista=Object.keys(D.viste)[0],sel=null,giro=null;
const cv=Q('c'),cx=cv.getContext('2d'),img=new Image(),off=document.createElement('canvas');
function lettura(o,k,mod){const u=mod?o.MU[k]:o.U[k],w=mod?o.MW[k]:o.W[k];return{kn:mod?o.MS[k]:o.S[k],dir:(Math.atan2(-u,-w)*180/Math.PI+360)%360}}
function valore(o,k){if(strato==='mod')return o.MS[k];if(strato==='raff')return o.G[k];if(strato==='diff'){const m=o.MS[k];return m>1.5?o.S[k]/m:1}return o.S[k]}
function colore(v){if(strato==='diff'){if(v<1)return[42,120,214,Math.min(210,(1-v)/0.6*210)];return[217,72,43,Math.min(210,(v-1)/0.6*210)]}
 const c=wcol(v);return[c[0],c[1],c[2],130]}
function sovrapposizione(o,v){const x0=v[0],y0=v[1],sw=v[2],sh=v[3],KO=Math.max(1,Math.min(4,Math.floor(800/sw))),w=sw*KO,hh=sh*KO;
 off.width=w;off.height=hh;const oc=off.getContext('2d'),im=oc.createImageData(w,hh),d=im.data;
 const F=new Float32Array(sw*sh);for(let j=0;j<sh;j++)for(let i=0;i<sw;i++)F[j*sw+i]=valore(o,(y0+j)*nx+x0+i);
 for(let y=0;y<hh;y++){const gy=(y+.5)/KO-.5,j0=Math.floor(gy),fy=gy-j0,ja=Math.max(0,j0),jb=Math.min(sh-1,j0+1);
  for(let x=0;x<w;x++){const gx=(x+.5)/KO-.5,i0=Math.floor(gx),fx=gx-i0,ia=Math.max(0,i0),ib=Math.min(sw-1,i0+1);
   const val=(F[ja*sw+ia]*(1-fx)+F[ja*sw+ib]*fx)*(1-fy)+(F[jb*sw+ia]*(1-fx)+F[jb*sw+ib]*fx)*fy,c=colore(val),p=4*(y*w+x);
   d[p]=c[0];d[p+1]=c[1];d[p+2]=c[2];d[p+3]=c[3]}}
 oc.putImageData(im,0,0);return off}
function freccia(X,Y,kn,dir,L,col,dpr){
 // come in Shift: freccia sottile che indica dove va il vento, colorata per intensita'
 const a=(dir+180)*Math.PI/180,dx=Math.sin(a),dy=-Math.cos(a),nx_=-dy,ny_=dx,hl=L*0.36,hw=L*0.15;
 const tx=X+dx*L/2,ty=Y+dy*L/2,sx=X-dx*L/2,sy=Y-dy*L/2,bx=tx-dx*hl,by=ty-dy*hl;
 cx.lineCap='round';cx.lineJoin='round';
 for(const [lw,c] of [[1.3*dpr+2.2*dpr,'rgba(255,255,255,.55)'],[1.3*dpr,col]]){cx.strokeStyle=c;cx.lineWidth=lw;cx.beginPath();cx.moveTo(sx,sy);cx.lineTo(bx,by);cx.stroke()}
 cx.beginPath();cx.moveTo(tx,ty);cx.lineTo(bx+nx_*hw,by+ny_*hw);cx.lineTo(bx-nx_*hw,by-ny_*hw);cx.closePath();
 cx.strokeStyle='rgba(255,255,255,.55)';cx.lineWidth=2*dpr;cx.stroke();cx.fillStyle=col;cx.fill()}
function barba(X,Y,kn,dir,L,dpr){
 const a=dir*Math.PI/180,sx=Math.sin(a),sy=-Math.cos(a),rx=-sy,ry=sx,r5=Math.round(kn/5)*5,asta=new Path2D(),pieno=new Path2D();
 if(r5<5){asta.arc(X,Y,L*0.17,0,7)}
 else{const Tx=X+sx*L,Ty=Y+sy*L,fl=L*0.46,pas=L*0.18,ca=Math.cos(0.35),sa=Math.sin(0.35);asta.moveTo(X,Y);asta.lineTo(Tx,Ty);
  let n50=Math.floor(r5/50),resto=r5-50*n50,n10=Math.floor(resto/10),n5=resto-10*n10>=5?1:0,d=0;
  for(let i=0;i<n50;i++){const bx=Tx-sx*d,by=Ty-sy*d;pieno.moveTo(bx,by);pieno.lineTo(bx+rx*fl,by+ry*fl);pieno.lineTo(bx-sx*pas*1.4,by-sy*pas*1.4);pieno.closePath();d+=pas*1.7}
  for(let i=0;i<n10;i++){const bx=Tx-sx*d,by=Ty-sy*d;asta.moveTo(bx,by);asta.lineTo(bx+(rx*ca+sx*sa)*fl,by+(ry*ca+sy*sa)*fl);d+=pas}
  if(n5){if(n50+n10===0)d=pas;const bx=Tx-sx*d,by=Ty-sy*d;asta.moveTo(bx,by);asta.lineTo(bx+(rx*ca+sx*sa)*fl*0.5,by+(ry*ca+sy*sa)*fl*0.5)}}
 cx.lineCap='round';cx.lineJoin='round';cx.strokeStyle='rgba(255,255,255,.7)';cx.lineWidth=3.2*dpr;cx.stroke(asta);cx.stroke(pieno);
 cx.strokeStyle='#1B1B1B';cx.lineWidth=1.4*dpr;cx.stroke(asta);cx.fillStyle='#1B1B1B';cx.fill(pieno)}
function vento(o,v,px,dpr,css){const x0=v[0],y0=v[1],sw=v[2],sh=v[3],pxc=css/sw,mod=strato==='mod';
 // griglia fissa sulla carta: stessi punti a ogni ora; il passo dipende solo dallo zoom
 const passi=[2,3,4,6,8,12,16,24,32];let st=passi.find(p=>p*pxc>=(barbe?30:19))||32;
 const L=barbe?Math.min(st*px*0.78,24*dpr):Math.min(st*px*0.74,15*dpr),r=Math.min(2,Math.floor(st/3));
 for(let cj=Math.floor(st/2);cj<ny;cj+=st){if(cj<y0||cj>=y0+sh)continue;
  for(let ci=Math.floor(st/2);ci<nx;ci+=st){if(ci<x0||ci>=x0+sw)continue;
   let su=0,sv=0,ss=0,n=0;
   for(let dj=-r;dj<=r;dj++)for(let di=-r;di<=r;di++){const a=ci+di,b=cj+dj;if(a<0||b<0||a>=nx||b>=ny)continue;const q=b*nx+a;
    su+=mod?o.MU[q]:o.U[q];sv+=mod?o.MW[q]:o.W[q];ss+=mod?o.MS[q]:o.S[q];n++}
   const kn=ss/n,dir=(Math.atan2(-su,-sv)*180/Math.PI+360)%360,X=(ci-x0+.5)*px,Y=(cj-y0+.5)*px;
   if(barbe)barba(X,Y,kn,dir,L,dpr);
   else{if(kn<0.8){cx.fillStyle='#1F2D3A';cx.beginPath();cx.arc(X,Y,1.2*dpr,0,7);cx.fill();continue}
    freccia(X,Y,kn,dir,L,(campitura||strato==='diff'||strato==='raff')?'#1F2D3A':rgb(wcol(kn)),dpr)}}}}
function etichette(v,px,W,H,dpr){const x0=v[0],y0=v[1],sw=v[2],presi=[];
 cx.textBaseline='middle';
 const libero=r=>presi.every(q=>r[2]<q[0]||r[0]>q[2]||r[3]<q[1]||r[1]>q[3]);
 const minuti=sw<nx*0.35;
 for(const l of D.luoghi){if(l[3]>=4&&!minuti)continue;const X=(l[1]*nx-x0)*px,Y=(l[2]*ny-y0)*px;if(X<6||Y<10||X>W-6||Y>H-10)continue;
  const picco=l[4]==='p',fs=(l[3]<=1?12.5:11.5)*dpr;cx.font=(l[3]<=2?'600 ':'500 ')+fs+'px system-ui,sans-serif';
  // la scritta prova a destra, a sinistra, sopra e sotto il punto; se non trova posto libero non si scrive
  const tw=cx.measureText(l[0]).width,g=7*dpr,v1=fs*1.15,pos=[[g,0],[-tw-g,0],[-tw/2,-v1],[-tw/2,v1],[-g,-v1],[-tw+g,-v1],[-g,v1],[-tw+g,v1]];let dx=null,dy=0,r=null;
  for(const p of pos){let tx=X+p[0];const ty=Y+p[1];
   if(p[1]!==0)tx=Math.max(2*dpr,Math.min(W-tw-2*dpr,tx));   // sopra e sotto la scritta puo' scorrere per restare dentro la carta
   if(tx<2*dpr||tx+tw>W-2*dpr||ty<fs*0.8||ty>H-fs*0.8)continue;
   const rr=[Math.min(X,tx)-3*dpr,Math.min(Y,ty)-fs*0.75,Math.max(X,tx+tw)+3*dpr,Math.max(Y,ty)+fs*0.75];if(libero(rr)){dx=tx-X;dy=p[1];r=rr;break}}
  if(dx===null)continue;presi.push(r);
  if(picco){cx.fillStyle='#5a3d1e';cx.beginPath();cx.moveTo(X,Y-5*dpr);cx.lineTo(X+4.5*dpr,Y+3*dpr);cx.lineTo(X-4.5*dpr,Y+3*dpr);cx.closePath();cx.fill()}
  else{cx.fillStyle='#c0392b';cx.beginPath();cx.arc(X,Y,2.6*dpr,0,7);cx.fill()}
  cx.lineWidth=3.2*dpr;cx.strokeStyle='rgba(255,253,245,.92)';cx.strokeText(l[0],X+dx,Y+dy);cx.fillStyle=picco?'#4a3420':'#14140f';cx.fillText(l[0],X+dx,Y+dy)}}
function fumetto(X,Y,t1,t2,W,dpr){
 cx.font='700 '+(16*dpr)+'px system-ui,sans-serif';const w1=cx.measureText(t1).width;cx.font='500 '+(13*dpr)+'px system-ui,sans-serif';const w2=cx.measureText(t2).width;
 const bw=Math.max(w1,w2)+18*dpr,bh=46*dpr;let bx=Math.max(4*dpr,Math.min(W-bw-4*dpr,X-bw/2)),by=Y-bh-16*dpr;if(by<4*dpr)by=Y+16*dpr;
 cx.fillStyle='rgba(255,253,245,.96)';cx.strokeStyle='#1F2D3A';cx.lineWidth=1.5*dpr;cx.beginPath();if(cx.roundRect)cx.roundRect(bx,by,bw,bh,8*dpr);else cx.rect(bx,by,bw,bh);cx.fill();cx.stroke();
 cx.fillStyle='#1F2D3A';cx.textBaseline='alphabetic';cx.font='700 '+(16*dpr)+'px system-ui,sans-serif';cx.fillText(t1,bx+9*dpr,by+20*dpr);
 cx.font='500 '+(13*dpr)+'px system-ui,sans-serif';cx.fillText(t2,bx+9*dpr,by+38*dpr)}
function disegna(){const v=D.viste[vista],x0=v[0],y0=v[1],sw=v[2],sh=v[3],o=ora(h);
 const css=cv.clientWidth,dpr=Math.min(3,window.devicePixelRatio||1),W=Math.round(css*dpr),H=Math.round(css*sh/sw*dpr);
 if(cv.width!==W||cv.height!==H){cv.width=W;cv.height=H}
 cx.imageSmoothingEnabled=true;cx.imageSmoothingQuality='high';cx.drawImage(img,x0*KB,y0*KB,sw*KB,sh*KB,0,0,W,H);
 if(campitura||strato==='diff'||strato==='raff')cx.drawImage(sovrapposizione(o,v),0,0,W,H);
 const px=W/sw;
 if(strato==='ril'){cx.fillStyle='rgba(31,45,58,.45)';const st=Math.max(1,Math.round(sw/70));
  for(let j=0;j<sh;j+=st)for(let i=0;i<sw;i+=st){if(o.R[(y0+j)*nx+x0+i]){cx.beginPath();cx.arc((i+.5)*px,(j+.5)*px,Math.max(0.9*dpr,px*0.16),0,7);cx.fill()}}}
 vento(o,v,px,dpr,css);
 etichette(v,px,W,H,dpr);
 if(sel!=null){const si=sel%nx-x0,sj=Math.floor(sel/nx)-y0;
  if(si>=0&&sj>=0&&si<sw&&sj<sh){const X=(si+.5)*px,Y=(sj+.5)*px,w=lettura(o,sel,strato==='mod');
   cx.lineWidth=4*dpr;cx.strokeStyle='#fffdf5';cx.beginPath();cx.arc(X,Y,8*dpr,0,7);cx.stroke();cx.lineWidth=2*dpr;cx.strokeStyle='#1F2D3A';cx.stroke();
   fumetto(X,Y,num(w.kn,1)+' kn · '+forza(w.kn)+(strato!=='mod'?' · raff. '+num(o.G[sel]):''),'da '+Math.round(w.dir)+'° '+rosa(w.dir)+' · '+nome(w.dir),W,dpr)}}
 testi()}
function testi(){const o=ORE[h],t=new Date(o.t+':00Z');
 Q('tl').textContent=new Intl.DateTimeFormat('it',{weekday:'short',day:'numeric',month:'numeric',hour:'2-digit',minute:'2-digit',timeZone:'Europe/Rome'}).format(t);
 const st=a=>a>=1?'rimescolata':a>=0.25?'poco stabile':'stabile';
 const reg=o.a2>=1?'il vento scavalca i rilievi':o.a2>=0.25?'il vento in parte scavalca, in parte aggira i rilievi':'il vento aggira i rilievi, ridossi lunghi';
 Q('reg').textContent=(o.q?'A 750 m: '+o.q[0]+' kn da '+o.q[1]+'° · ':'')+'aria '+st(o.a2)+' fino alla cima dei rilievi'+(o.a22!==undefined&&o.a22!==o.a2?', '+st(o.a22)+' sopra':'')+': '+reg;
 if(strato==='diff')Q('leg').innerHTML='<div class="scala" style="background:linear-gradient(90deg,rgb(42,120,214),rgba(42,120,214,0) 50%,rgba(217,72,43,0) 50%,rgb(217,72,43))"></div><div class="tacche"><span>−60%</span><span>ridosso</span><span>come il modello</span><span>accelerazione</span><span>+60%</span></div>';
 else Q('leg').innerHTML='<div class="scala" style="background:linear-gradient(90deg,'+WS.map(s=>rgb(s[1])+' '+(s[0]/32*100).toFixed(0)+'%').join(',')+')"></div><div class="tacche">'+WS.map(s=>'<span>'+s[0]+'</span>').join('')+'</div><div class="reg">'+(strato==='raff'?'raffiche, ':'')+'nodi a 10 m dalla superficie'+(strato==='ril'?' · puntini: ricircolo probabile':'')+'</div>';
 if(sel==null){Q('info').textContent='Tocca la carta: nel punto compaiono intensità e direzione.';return}
 const o2=ora(h),r=lettura(o2,sel,false),m=lettura(o2,sel,true),q=m.kn>1.5?(r.kn/m.kn-1)*100:0;
 Q('info').innerHTML='<b>'+num(r.kn,1)+' kn da '+Math.round(r.dir)+'° ('+rosa(r.dir)+', '+nome(r.dir)+')</b>, raffiche '+num(o2.G[sel])+' kn'+(o2.R[sel]?' · <b>ricircolo probabile</b>':'')+
  '<br>modello a 2,2 km: '+num(m.kn,1)+' kn da '+Math.round(m.dir)+'° · differenza '+(q>=0?'+':'')+num(q)+'% · '+(MAR[sel]?'mare':'terra, quota '+quota[sel]+' m')}
function premi(gruppo,valore){document.querySelectorAll('#'+gruppo+' button[data-v]').forEach(b=>b.setAttribute('aria-pressed',b.dataset.v===valore))}
Q('viste').innerHTML=Object.keys(D.viste).map(v=>'<button type="button" data-v="'+v+'">'+v+'</button>').join('');
Q('viste').addEventListener('click',e=>{const b=e.target.closest('button');if(!b)return;vista=b.dataset.v;premi('viste',vista);disegna()});
Q('strati').addEventListener('click',e=>{const b=e.target.closest('button');if(!b)return;
 if(b.id==='camp'){campitura=!campitura;b.setAttribute('aria-pressed',campitura)}
 else if(b.id==='barbe'){barbe=!barbe;b.setAttribute('aria-pressed',barbe)}
 else{strato=b.dataset.v;premi('strati',strato)}disegna()});
Q('cur').max=ORE.length-1;Q('cur').addEventListener('input',e=>{h=+e.target.value;disegna()});
function vai(d){h=Math.max(0,Math.min(ORE.length-1,h+d));Q('cur').value=h;disegna()}
Q('prima').addEventListener('click',()=>vai(-1));Q('dopo').addEventListener('click',()=>vai(1));
Q('play').addEventListener('click',()=>{if(giro){clearInterval(giro);giro=null;Q('play').textContent='▶';return}
 Q('play').textContent='❚❚';giro=setInterval(()=>{h=(h+1)%ORE.length;Q('cur').value=h;disegna()},900)});
cv.addEventListener('click',e=>{const r=cv.getBoundingClientRect(),v=D.viste[vista],i=Math.min(v[2]-1,Math.max(0,Math.floor((e.clientX-r.left)/r.width*v[2]))),j=Math.min(v[3]-1,Math.max(0,Math.floor((e.clientY-r.top)/r.height*v[3])));sel=(v[1]+j)*nx+v[0]+i;disegna()});
window.addEventListener('resize',disegna);
premi('viste',vista);premi('strati',strato);
img.onload=disegna;img.src=Q('base').src;
})();
"""


ESEMPIO = """<svg width="150" height="34" viewBox="0 0 150 34" aria-hidden="true"><g fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round">
<path d="M8 28H62M62 28l5 -14M55 28l2.5 -7"/><circle cx="8" cy="28" r="1.8" fill="currentColor"/></g>
<text x="72" y="31" font-size="11" fill="currentColor">15 nodi, da ovest</text></svg>"""


def pagina(dati, base, mime="image/webp"):
    import datetime as dt
    from zoneinfo import ZoneInfo
    gen = dt.datetime.fromisoformat(dati["generato"]).astimezone(ZoneInfo("Europe/Rome"))
    modello = dati.get("nome_modello") or dati["modello"]
    titolo = dati.get("titolo") or "Vento ad alta definizione"
    return f"""<!DOCTYPE html>
<html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="robots" content="noindex,nofollow"><title>{titolo}</title><style>{CSS}</style></head><body><main>
<h1>{titolo}</h1>
<p class="sub">Vento a 10 m · dettaglio di {dati['dx']:.0f} m · {modello} riportato sul rilievo · {dati['n_ore']} ore · calcolato {gen:%d/%m alle %H:%M}</p>
<div class="card">
<div class="riga" id="viste"></div>
<div class="riga" id="strati"><button type="button" data-v="ril">Sul rilievo</button><button type="button" data-v="mod">Modello a 2,2 km</button><button type="button" data-v="raff">Raffiche</button><button type="button" data-v="diff">Differenza</button><button type="button" id="camp" aria-pressed="false">Campitura</button><button type="button" id="barbe" aria-pressed="false">Banderuole</button></div>
<canvas id="c" width="600" height="600"></canvas>
<div class="tempo"><button type="button" id="prima" aria-label="ora prima">◀</button><button type="button" id="play" aria-label="scorri">▶</button><button type="button" id="dopo" aria-label="ora dopo">▶|</button><input type="range" id="cur" min="0" max="14" value="0" aria-label="ora"></div>
<div class="tl" id="tl"></div><div class="reg" id="reg"></div>
<div id="leg"></div><p class="info" id="info"></p>
</div>
<div class="card"><ul class="note">
<li><b>Frecce</b>: indicano dove va il vento, colorate per intensità con la scala qui sopra. Stanno sempre negli stessi punti: scorrendo le ore girano e cambiano colore. <b>Banderuole</b> le sostituisce con le barbe: l'asta punta da dove viene il vento, tratto lungo 10 nodi, corto 5. {ESEMPIO}</li>
<li><b>Campitura</b>: colora tutta la carta con l'intensità del vento, per vedere a colpo d'occhio ridossi e accelerazioni. <b>Raffiche</b>: la raffica attesa in ogni punto, più forte dove il vento arriva da terra o esce dalla scia di un rilievo; è tarata sul rapporto tra raffica e vento medio che il modello dà al largo.</li>
<li><b>Dettaglio di {dati['dx']:.0f} m</b>: la carta è divisa in quadrati di {dati['dx']:.0f} m di lato, ognuno con il suo vento. Non è un'altezza: il vento è sempre a 10 m dalla superficie, cioè a 10 m sul mare in mare, come nei GRIB e negli anemometri. Il calcolo comprende un margine attorno alla carta e usa l'aria fino a 3000 m, perché decide come il vento passa sopra o attorno ai rilievi.</li>
<li><b>Sul rilievo</b>: il vento del modello, con il suo profilo verticale fino a 3000 m, ricalcolato in tre dimensioni perché non attraversi le montagne e conservi la massa, con la stabilità dell'aria presa dal modello in due strati (sotto e sopra la cima dei rilievi), più due stime: il ridosso sottovento ai rilievi e la frenata del vento che arriva da terra, che sul mare riaccelera in una fascia lungo costa. <b>Differenza</b>: dove il rilievo toglie vento rispetto al modello (blu) e dove lo aggiunge (rosso).</li>
<li><b>Puntini</b>: rilievo ripido subito sopravvento, ricircolo probabile. Indicano la zona, non la forma del vortice.</li>
<li>Lunghezza e forza dei ridossi e della fascia lungo costa sono stime non ancora confrontate con anemometri veri. Brezze di versante e termiche non sono calcolate: c'è solo la brezza che dà il modello.</li>
</ul><p class="note">Fonti: {modello} e profilo verticale via Open-Meteo.com, CC BY 4.0; rilievo Copernicus DEM GLO-30; copertura del suolo ESA WorldCover 2021, CC BY 4.0; fiumi e nomi © OpenStreetMap, ODbL.</p></div>
<img id="base" alt="" hidden src="data:{mime};base64,{base}">
</main><script>{JS.replace("__DATI__", json.dumps(dati, separators=(",", ":"), ensure_ascii=False))}</script></body></html>"""
