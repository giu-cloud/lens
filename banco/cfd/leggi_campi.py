import numpy as np
def leggi(p, n_venti):
    b=open(p,'rb').read()
    nx,ny=np.frombuffer(b,'>i4',2,0); box=np.frombuffer(b,'>f8',4,8)
    o=8+32; N=nx*ny; out=[]
    for _ in range(n_venti):
        coppia=[]
        for _ in range(2):
            u=np.frombuffer(b,'>f4',N,o).reshape(ny,nx); o+=4*N; v=np.frombuffer(b,'>f4',N,o).reshape(ny,nx); o+=4*N
            coppia.append((u,v))
        out.append(coppia)
    s,n,w,e=box
    lat=n-(np.arange(ny)+0.5)*(n-s)/ny; lon=w+(np.arange(nx)+0.5)*(e-w)/nx
    return lat,lon,out
