"""Un caso di fluidodinamica vera (OpenFOAM, RANS k-epsilon, aria neutra) sul rilievo di un atlante di Lens, per
confrontarlo con Lens (conservazione della massa) su una direzione del vento.

  python3 caso_cfd.py <atlante .lna> <cartella del caso> --centro LAT,LON --da GRADI [--vento M/S]
         [--x -14,10] [--y -10,8] [--dx 150] [--alto 3000] [--strati 30] [--primo 10]

Il dominio è ruotato: l'asse x va dove va il vento (ingresso a ovest per un vento da ovest), y a sinistra. Il rilievo e
la scabrezza vengono dall'atlante (Copernicus GLO-30 e WorldCover, già ricampionati da Lens); ai bordi il rilievo
sfuma a zero su 2 km, così l'ingresso e l'uscita stanno su terreno piatto. Maglia: blockMesh (un blocco), poi i punti
si alzano sul rilievo in coordinate che seguono il terreno. Al suolo la funzione di parete per terreno scabro con la
z0 di ogni faccia; all'ingresso il profilo dello strato limite neutro (Richards e Hoxey).
Dopo: blockMesh, python3 caso_cfd.py ... --solo-punti, checkMesh, simpleFoam (vedi esegui.sh).
"""
import argparse
import gzip
import math
import os
import struct
import sys

import numpy as np


def leggi_atlante(p):
    b = gzip.open(p).read(); o = 4
    NX, NY, nx, ny, m, k = struct.unpack_from("<6i", b, o); o += 24
    S, N, O, E, s, n, w, e = struct.unpack_from("<8d", b, o); o += 64
    dx, H, Neff = struct.unpack_from("<3f", b, o); o += 12; o += 4 * k
    h = np.frombuffer(b, "<i2", NX * NY, o).reshape(NY, NX).astype("f8"); o += 2 * NX * NY
    z0 = np.frombuffer(b, "<f4", NX * NY, o).reshape(NY, NX).astype("f8"); o += 4 * NX * NY
    mare = np.frombuffer(b, "u1", NX * NY, o).reshape(NY, NX).astype(bool)
    # righe da nord: lat della riga r = N - (r + 0,5) * passo
    lat = N - (np.arange(NY) + 0.5) * (N - S) / NY
    lon = O + (np.arange(NX) + 0.5) * (E - O) / NX
    return dict(h=h, z0=z0, mare=mare, lat=lat, lon=lon)


def campiona(A, la, lo, campo):
    """Bilineare sulla griglia dell'atlante (righe da nord)."""
    lat, lon = A["lat"], A["lon"]
    fy = (lat[0] - la) / (lat[0] - lat[1]); fx = (lo - lon[0]) / (lon[1] - lon[0])
    j0 = np.clip(np.floor(fy).astype(int), 0, len(lat) - 2); i0 = np.clip(np.floor(fx).astype(int), 0, len(lon) - 2)
    ty = np.clip(fy - j0, 0, 1); tx = np.clip(fx - i0, 0, 1)
    F = A[campo]
    return ((1 - ty) * (1 - tx) * F[j0, i0] + (1 - ty) * tx * F[j0, i0 + 1] + ty * (1 - tx) * F[j0 + 1, i0] + ty * tx * F[j0 + 1, i0 + 1])


class Geometria:
    def __init__(self, lat0, lon0, da_gradi):
        self.lat0, self.lon0 = lat0, lon0
        verso = math.radians(da_gradi + 180.0)            # dove va il vento, dal nord in senso orario
        self.f = np.array([math.sin(verso), math.cos(verso)])   # (est, nord)
        self.l = np.array([-self.f[1], self.f[0]])              # a sinistra del vento
        self.kx = 111.19 * math.cos(math.radians(lat0)); self.ky = 111.19

    def latlon(self, xr, yr):
        """km nel dominio ruotato -> lat, lon."""
        e = xr * self.f[0] + yr * self.l[0]; n = xr * self.f[1] + yr * self.l[1]
        return self.lat0 + n / self.ky, self.lon0 + e / self.kx

    def ruotato(self, la, lo):
        e = (lo - self.lon0) * self.kx; n = (la - self.lat0) * self.ky
        return e * self.f[0] + n * self.f[1], e * self.l[0] + n * self.l[1]


def griglia_z(alto, strati, primo):
    """Spessori geometrici: il primo 'primo' m, la somma 'alto'. Restituisce il rapporto tra ultimo e primo."""
    lo, hi = 1.0001, 2.0
    for _ in range(100):
        q = 0.5 * (lo + hi)
        if primo * (q ** strati - 1) / (q - 1) > alto:
            hi = q
        else:
            lo = q
    return q ** (strati - 1)


def scrivi(p, testo):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(testo)


INTESTA = """FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}
"""


def caso(a, A, G):
    x0, x1 = a.x; y0, y1 = a.y
    nx = int(round((x1 - x0) * 1000 / a.dx)); ny = int(round((y1 - y0) * 1000 / a.dx))
    grad = griglia_z(a.alto, a.strati, a.primo)
    c = a.caso
    scrivi(f"{c}/system/blockMeshDict", INTESTA.format(cls="dictionary", obj="blockMeshDict") + f"""
convertToMeters 1;
vertices ( ({x0*1000} {y0*1000} 0) ({x1*1000} {y0*1000} 0) ({x1*1000} {y1*1000} 0) ({x0*1000} {y1*1000} 0)
           ({x0*1000} {y0*1000} {a.alto}) ({x1*1000} {y0*1000} {a.alto}) ({x1*1000} {y1*1000} {a.alto}) ({x0*1000} {y1*1000} {a.alto}) );
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {a.strati}) simpleGrading (1 1 {grad:.6f}) );
boundary (
  inlet  {{ type patch; faces ((0 4 7 3)); }}
  outlet {{ type patch; faces ((1 2 6 5)); }}
  lati   {{ type patch; faces ((0 1 5 4) (3 7 6 2)); }}
  top    {{ type patch; faces ((4 5 6 7)); }}
  ground {{ type wall;  faces ((0 3 2 1)); }}
);
""")
    scrivi(f"{c}/system/controlDict", INTESTA.format(cls="dictionary", obj="controlDict") + f"""
application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime {a.iter}; deltaT 1;
writeControl timeStep; writeInterval {a.iter}; purgeWrite 1; writeFormat ascii; writePrecision 6; writeCompression off;
timeFormat general; runTimeModifiable true; libs ("libatmosphericModels.so");
""")
    scrivi(f"{c}/system/fvSchemes", INTESTA.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default steadyState; }
gradSchemes { default cellLimited Gauss linear 1; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,k) bounded Gauss upwind;
  div(phi,epsilon) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }
wallDist { method meshWave; }
""")
    scrivi(f"{c}/system/fvSolution", INTESTA.format(cls="dictionary", obj="fvSolution") + """
solvers {
  p { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.05; }
  "(U|k|epsilon)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }
}
SIMPLE { nNonOrthogonalCorrectors 1; consistent yes; residualControl { p 2e-4; U 2e-5; "(k|epsilon)" 2e-5; } }
relaxationFactors { equations { U 0.7; ".*" 0.5; } fields { p 0.7; } }
""")
    scrivi(f"{c}/system/decomposeParDict", INTESTA.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {a.processori}; method simple; coeffs {{ n ({a.processori} 1 1); }}
""")
    scrivi(f"{c}/constant/transportProperties", INTESTA.format(cls="dictionary", obj="transportProperties") + "transportModel Newtonian; nu 1.5e-05;\n")
    scrivi(f"{c}/constant/turbulenceProperties", INTESTA.format(cls="dictionary", obj="turbulenceProperties") + """
simulationType RAS;
RAS { RASModel kEpsilon; turbulence on; printCoeffs on;
  kEpsilonCoeffs { Cmu 0.09; C1 1.44; C2 1.92; C3 0; sigmak 1.0; sigmaEps 1.11; } }
""")
    return nx, ny


def campi_iniziali(a, z0_suolo):
    c = a.caso
    abl = f"flowDir (1 0 0); zDir (0 0 1); Uref {a.vento}; Zref 10; z0 uniform {a.z0_ingresso}; zGround uniform 0; value $internalField;"
    z0l = "nonuniform List<scalar> " + str(len(z0_suolo)) + "(" + " ".join(f"{v:.5g}" for v in z0_suolo) + ")"
    scrivi(f"{c}/0/U", INTESTA.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0]; internalField uniform ({a.vento} 0 0);
boundaryField {{
  inlet {{ type atmBoundaryLayerInletVelocity; {abl} }}
  outlet {{ type inletOutlet; inletValue uniform (0 0 0); value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type noSlip; }}
}}
""")
    k0 = a.vento ** 2 * 0.0025
    scrivi(f"{c}/0/k", INTESTA.format(cls="volScalarField", obj="k") + f"""
dimensions [0 2 -2 0 0 0 0]; internalField uniform {k0:.4g};
boundaryField {{
  inlet {{ type atmBoundaryLayerInletK; {abl} }}
  outlet {{ type inletOutlet; inletValue $internalField; value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type kqRWallFunction; value $internalField; }}
}}
""")
    scrivi(f"{c}/0/epsilon", INTESTA.format(cls="volScalarField", obj="epsilon") + f"""
dimensions [0 2 -3 0 0 0 0]; internalField uniform 0.01;
boundaryField {{
  inlet {{ type atmBoundaryLayerInletEpsilon; {abl} }}
  outlet {{ type inletOutlet; inletValue $internalField; value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type epsilonWallFunction; Cmu 0.09; kappa 0.41; E 9.8; value $internalField; }}
}}
""")
    scrivi(f"{c}/0/nut", INTESTA.format(cls="volScalarField", obj="nut") + f"""
dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;
boundaryField {{
  inlet {{ type calculated; value uniform 0; }}
  outlet {{ type calculated; value uniform 0; }}
  lati {{ type calculated; value uniform 0; }} top {{ type calculated; value uniform 0; }}
  ground {{ type nutkAtmRoughWallFunction; Cmu 0.09; kappa 0.41; E 9.8; z0 {z0l}; value uniform 0; }}
}}
""")
    scrivi(f"{c}/0/p", INTESTA.format(cls="volScalarField", obj="p") + """
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;
boundaryField { inlet { type zeroGradient; } outlet { type fixedValue; value uniform 0; }
  lati { type slip; } top { type slip; } ground { type zeroGradient; } }
""")


def leggi_lista(p, tipo=float):
    """Una lista di OpenFOAM in ASCII: 'n' a capo '(' ... ')'. Restituisce n e il testo tra le parentesi."""
    import re
    t = open(p).read()
    m = re.search(r"\n(\d+)\s*\n\(", t)
    n = int(m.group(1))
    corpo = t[m.end():t.rindex(")")]
    return n, corpo


def alza_punti(a, A, G):
    """Riscrive constant/polyMesh/points sul rilievo e restituisce la z0 delle facce del suolo."""
    c = a.caso
    pp = f"{c}/constant/polyMesh/points"
    n, corpo = leggi_lista(pp)
    P = np.array(corpo.replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)
    xr, yr = P[:, 0] / 1000.0, P[:, 1] / 1000.0
    la, lo = G.latlon(xr, yr)
    h = np.maximum(campiona(A, la, lo, "h"), 0.0)
    # sfumatura ai bordi del dominio, su 2 km
    x0, x1 = a.x; y0, y1 = a.y
    d = np.minimum.reduce([xr - x0, x1 - xr, yr - y0, y1 - yr])
    w = np.clip(d / 2.0, 0, 1); w = 0.5 - 0.5 * np.cos(np.pi * w)
    h = h * w
    z = P[:, 2]
    P[:, 2] = h + z * (a.alto - h) / a.alto
    with open(pp, "w") as f:
        f.write(INTESTA.format(cls="vectorField", obj="points") + f"\n{n}\n(\n")
        for x, y, zz in P:
            f.write(f"({x:.3f} {y:.3f} {zz:.3f})\n")
        f.write(")\n")
    # facce del suolo: dal file boundary, poi i loro centri
    b = open(f"{c}/constant/polyMesh/boundary").read()
    k = b.index("ground"); seg = b[k:b.index("}", k)]
    nf = int(seg.split("nFaces")[1].split(";")[0]); st = int(seg.split("startFace")[1].split(";")[0])
    nfaces, corpoF = leggi_lista(f"{c}/constant/polyMesh/faces")
    righe = [r for r in corpoF.split("\n") if "(" in r]
    cen = []
    for r in righe[st:st + nf]:
        ids = [int(x) for x in r[r.index("(") + 1:r.index(")")].split()]
        cen.append(P[ids].mean(axis=0))
    cen = np.array(cen)
    la, lo = G.latlon(cen[:, 0] / 1000.0, cen[:, 1] / 1000.0)
    z0 = np.clip(campiona(A, la, lo, "z0"), 0.0002, 1.0)
    mare = campiona(A, la, lo, "mare") > 0.5
    z0[mare] = 0.0002
    return z0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlante"); ap.add_argument("caso")
    ap.add_argument("--centro", required=True); ap.add_argument("--da", type=float, required=True)
    ap.add_argument("--vento", type=float, default=8.0); ap.add_argument("--z0-ingresso", type=float, default=0.0002)
    ap.add_argument("--x", default="-14,10"); ap.add_argument("--y", default="-10,8")
    ap.add_argument("--dx", type=float, default=150.0); ap.add_argument("--alto", type=float, default=3000.0)
    ap.add_argument("--strati", type=int, default=30); ap.add_argument("--primo", type=float, default=10.0)
    ap.add_argument("--iter", type=int, default=1500); ap.add_argument("--processori", type=int, default=2)
    ap.add_argument("--solo-punti", action="store_true")
    a = ap.parse_args()
    a.x = [float(v) for v in a.x.split(",")]; a.y = [float(v) for v in a.y.split(",")]
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    A = leggi_atlante(a.atlante); G = Geometria(lat0, lon0, a.da)
    if not a.solo_punti:
        nx, ny = caso(a, A, G)
        print(f"caso scritto: {nx} x {ny} x {a.strati} celle; ora blockMesh, poi --solo-punti")
        return
    z0 = alza_punti(a, A, G)
    campi_iniziali(a, z0)
    print(f"punti sul rilievo, {len(z0)} facce di suolo (z0 mediana {np.median(z0):.4f})")


if __name__ == "__main__":
    main()
