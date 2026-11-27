"""CFD con l'aria stratificata (27/11): come caso_cfd.py, ma con il solutore di Boussinesq (buoyantBoussinesqSimpleFoam)
e la temperatura potenziale. Piano e criterio in lens/cfd_piano_26_11.md.

  python3 caso_cfd_s.py <terreno> <caso> --centro LAT,LON --da GRADI --N 0.01 [--vento 8] [--x -22,22] [--y -22,22]
         [--dx 400] [--alto 10000] [--strati 40] [--primo 10] [--delta 500] [--spugna 6000] [--iter 2000] [--processori 4]
  passi: caso_cfd_s.py ...; blockMesh; caso_cfd_s.py ... --solo-punti; checkMesh; topoSet;
         caso_cfd_s.py ... --stratifica; buoyantBoussinesqSimpleFoam

Che cosa cambia rispetto alla CFD neutra del 25/11:
- l'aria ha N costante (frequenza di Brunt-Vaisala): temperatura potenziale theta = TRef + Gamma z, Gamma = N^2 TRef / g.
  N = 0 è la stessa catena in aria neutra (il controllo G-N);
- all'ingresso uno strato limite alto delta: vento logaritmico fino a delta e poi costante, turbolenza che si spegne a
  delta (la turbolenza del profilo di Richards e Hoxey, costante fino in cima, mescolerebbe tutta la colonna e
  cancellerebbe la stratificazione);
- il tetto è alto (10 km) e sopra 'spugna' c'è uno strato che riporta vento e temperatura a quelli d'ingresso
  (smorzamento di Rayleigh, fvOptions semiImplicitSource), perché le onde di gravità non rimbalzino;
- all'uscita la pressione è quella idrostatica dell'aria stratificata;
- il rilievo è lisciato alla scala della maglia (gaussiana, sigma = dx) e sfumato su 6 km verso i bordi: con i pendii
  ripidi di Peloritani e Aspromonte e le celle basse (10 m) il solutore di Boussinesq di OpenFOAM 1912 esplode nelle
  prime iterazioni (simpleFoam no); con il rilievo alla scala della maglia converge;
- la convergenza si giudica sul campo a 10 m a iter-300 e a iter (il flusso cfd_s.yml fa due tratti).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from caso_cfd import Geometria, INTESTA, campiona, griglia_z, leggi_atlante, leggi_lista, scrivi   # noqa: E402

G0, TREF, KAPPA, CMU = 9.81, 290.0, 0.41, 0.09
BETA = 1.0 / TREF
ZONE_SPUGNA = 4
# punti di controllo per sito: nome, lat, lon (vento a 10 m)
PUNTI = {
    "stretto": [("imbocco_nord", 38.29, 15.67), ("ganzirri", 38.245, 15.63), ("messina_villa", 38.20, 15.615),
                ("messina_reggio", 38.12, 15.59), ("sud_reggio", 38.09, 15.585), ("tirreno_largo", 38.40, 15.55),
                ("ionio_largo", 37.98, 15.45)],
}


def alza_punti(a, A, G):
    """Come caso_cfd.alza_punti, ma la sfumatura del rilievo verso i bordi è larga (--sfuma km, 6 di solito): con 2 km
    il rilievo tagliato al bordo diventa un dirupo di 40 gradi e il solutore di Boussinesq (che ricostruisce il gradiente
    di pressione dalle facce) esplode lì nelle prime iterazioni. Restituisce la z0 delle facce del suolo."""
    c = a.caso
    pp = f"{c}/constant/polyMesh/points"
    n, corpo = leggi_lista(pp)
    P = np.array(corpo.replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)
    xr, yr = P[:, 0] / 1000.0, P[:, 1] / 1000.0
    la, lo = G.latlon(xr, yr)
    if a.liscia > 0:            # il rilievo lisciato alla scala della maglia (gaussiana, sigma in m)
        from scipy.ndimage import gaussian_filter
        passo = abs(A["lat"][0] - A["lat"][1]) * 111190.0
        A = dict(A, h=gaussian_filter(A["h"], a.liscia / passo))
    h = np.maximum(campiona(A, la, lo, "h"), 0.0)
    x0, x1 = a.x; y0, y1 = a.y
    d = np.minimum.reduce([xr - x0, x1 - xr, yr - y0, y1 - yr])
    w = np.clip(d / a.sfuma, 0, 1); w = 0.5 - 0.5 * np.cos(np.pi * w)
    h = h * w
    z = P[:, 2]
    P[:, 2] = h + z * (a.alto - h) / a.alto
    with open(pp, "w") as f:
        f.write(INTESTA.format(cls="vectorField", obj="points") + f"\n{n}\n(\n")
        for x, y, zz in P:
            f.write(f"({x:.3f} {y:.3f} {zz:.3f})\n")
        f.write(")\n")
    b = open(f"{c}/constant/polyMesh/boundary").read()
    k = b.index("ground"); seg = b[k:b.index("}", k)]
    nf = int(seg.split("nFaces")[1].split(";")[0]); st = int(seg.split("startFace")[1].split(";")[0])
    nfaces, corpoF = leggi_lista(f"{c}/constant/polyMesh/faces")
    righe = [r for r in corpoF.split("\n") if "(" in r]
    cen = np.array([P[[int(x) for x in r[r.index("(") + 1:r.index(")")].split()]].mean(axis=0) for r in righe[st:st + nf]])
    la, lo = G.latlon(cen[:, 0] / 1000.0, cen[:, 1] / 1000.0)
    z0 = np.clip(campiona(A, la, lo, "z0"), 0.0002, 1.0)
    z0[campiona(A, la, lo, "mare") > 0.5] = 0.0002
    return z0


def gamma(N):
    return N * N * TREF / G0


def profili(a):
    """Profili d'ingresso: z (m), vento, k, epsilon, theta, p_rgh idrostatica (cinematica)."""
    z0 = a.z0_ingresso
    ust = a.vento * KAPPA / math.log((10.0 + z0) / z0)
    z = np.unique(np.r_[0.0, np.geomspace(0.5, a.alto, 160), a.alto])
    zz = np.minimum(z, a.delta)
    U = ust / KAPPA * np.log((zz + z0) / z0)
    f = np.clip(1.0 - z / a.delta, 0.0, 1.0)
    k = ust ** 2 / math.sqrt(CMU) * f ** 2 + 1e-4
    lmix = KAPPA * (z + z0) / (1.0 + KAPPA * (z + z0) / 40.0)     # lunghezza di mescolamento con tetto di 40 m (Blackadar)
    eps = CMU ** 0.75 * k ** 1.5 / lmix
    Gm = gamma(a.N)
    th = TREF + Gm * z
    prgh = -G0 * BETA * Gm * z ** 2 / 2.0
    return z, U, k, eps, th, prgh, ust


def tabella(z, v, vettore=False):
    if vettore:
        return "table (" + " ".join(f"({zz:.3f} ({vv:.5g} 0 0))" for zz, vv in zip(z, v)) + ")"
    return "table (" + " ".join(f"({zz:.3f} {vv:.7g})" for zz, vv in zip(z, v)) + ")"


def caso(a, G):
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
application buoyantBoussinesqSimpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime {a.iter}; deltaT 1;
writeControl timeStep; writeInterval {a.iter}; purgeWrite 1; writeFormat ascii; writePrecision 6; writeCompression off;
timeFormat general; runTimeModifiable true; libs ("libatmosphericModels.so");
""")
    scrivi(f"{c}/system/fvSchemes", INTESTA.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default steadyState; }
gradSchemes { default cellLimited Gauss linear 1; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,T) bounded Gauss linearUpwind grad(T);
  div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }
wallDist { method meshWave; }
""")
    scrivi(f"{c}/system/fvSolution", INTESTA.format(cls="dictionary", obj="fvSolution") + """
solvers {
  p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.05; }
  "(U|T|k|epsilon)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
}
SIMPLE { momentumPredictor yes; nNonOrthogonalCorrectors 1; residualControl { p_rgh 1e-5; U 1e-6; T 1e-7; "(k|epsilon)" 1e-6; } }
relaxationFactors { fields { p_rgh 0.3; } equations { U 0.7; T 0.5; "(k|epsilon)" 0.5; } }
""")
    scrivi(f"{c}/system/decomposeParDict", INTESTA.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {a.processori}; method simple; coeffs {{ n ({a.processori} 1 1); }}
""")
    # zone della spugna, per quota (in alto il rilievo è quasi sparito: le coordinate seguono il terreno)
    zs = np.linspace(a.spugna, a.alto, ZONE_SPUGNA + 1)
    azioni = []
    for k in range(ZONE_SPUGNA):
        azioni.append(f"{{ name spugna{k}; type cellSet; action new; source boxToCell; sourceInfo {{ box (-1e9 -1e9 {zs[k]:.1f}) (1e9 1e9 {zs[k + 1]:.1f}); }} }}")
        azioni.append(f"{{ name spugna{k}; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set spugna{k}; }} }}")
    scrivi(f"{c}/system/topoSetDict", INTESTA.format(cls="dictionary", obj="topoSetDict") + "actions (\n" + "\n".join(azioni) + "\n);\n")
    scrivi(f"{c}/constant/g", INTESTA.format(cls="uniformDimensionedVectorField", obj="g") + "dimensions [0 1 -2 0 0 0 0]; value (0 0 -9.81);\n")
    scrivi(f"{c}/constant/transportProperties", INTESTA.format(cls="dictionary", obj="transportProperties") +
           f"transportModel Newtonian; nu 1.5e-05; beta {BETA:.6g}; TRef {TREF}; Pr 0.71; Prt 0.85;\n")
    scrivi(f"{c}/constant/turbulenceProperties", INTESTA.format(cls="dictionary", obj="turbulenceProperties") + """
simulationType RAS;
RAS { RASModel kEpsilon; turbulence on; printCoeffs on;
  kEpsilonCoeffs { Cmu 0.09; C1 1.44; C2 1.92; C3 0; sigmak 1.0; sigmaEps 1.11; } }
""")
    # la spugna: S = alfa (valore d'ingresso - campo), alfa da 0,002 a 0,02 1/s salendo
    z, U, k_, e_, th, prgh, ust = profili(a)
    Utop = float(U[-1])
    src = []
    for k in range(ZONE_SPUGNA):
        al = 0.002 * (10 ** (k / max(ZONE_SPUGNA - 1, 1)))
        zm = 0.5 * (zs[k] + zs[k + 1]); thm = TREF + gamma(a.N) * zm
        src.append(f"""spugnaU{k} {{ type vectorSemiImplicitSource; active true; selectionMode cellZone; cellZone spugna{k};
  volumeMode specific; injectionRateSuSp {{ U (({al * Utop:.6g} 0 0) {-al:.6g}); }} }}""")
        src.append(f"""spugnaT{k} {{ type scalarSemiImplicitSource; active true; selectionMode cellZone; cellZone spugna{k};
  volumeMode specific; injectionRateSuSp {{ T ({al * thm:.6g} {-al:.6g}); }} }}""")
    scrivi(f"{c}/constant/fvOptions", INTESTA.format(cls="dictionary", obj="fvOptions") + "\n".join(src) + "\n")
    return nx, ny


def campi(a, z0_suolo):
    c = a.caso
    z, U, k, eps, th, prgh, ust = profili(a)
    z0l = "nonuniform List<scalar> " + str(len(z0_suolo)) + "(" + " ".join(f"{v:.5g}" for v in z0_suolo) + ")"
    prof = lambda v, vet=False: f"type fixedProfile; profile {tabella(z, v, vet)}; direction (0 0 1); origin 0;"
    Gm = gamma(a.N)
    scrivi(f"{c}/0/U", INTESTA.format(cls="volVectorField", obj="U") + f"""
// Uref {a.vento} (vento a 10 m sul mare all'ingresso; campo_cfd.py lo legge da qui)
dimensions [0 1 -1 0 0 0 0]; internalField uniform ({float(U[-1]) * 0.8:.4g} 0 0);
boundaryField {{
  inlet {{ {prof(U, True)} }}
  outlet {{ type inletOutlet; inletValue uniform (0 0 0); value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type noSlip; }}
}}
""")
    scrivi(f"{c}/0/k", INTESTA.format(cls="volScalarField", obj="k") + f"""
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0.01;
boundaryField {{
  inlet {{ {prof(k)} }}
  outlet {{ type inletOutlet; inletValue uniform 1e-4; value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type kqRWallFunction; value $internalField; }}
}}
""")
    scrivi(f"{c}/0/epsilon", INTESTA.format(cls="volScalarField", obj="epsilon") + f"""
dimensions [0 2 -3 0 0 0 0]; internalField uniform 1e-4;
boundaryField {{
  inlet {{ {prof(eps)} }}
  outlet {{ type inletOutlet; inletValue uniform 1e-6; value $internalField; }}
  lati {{ type slip; }} top {{ type slip; }}
  ground {{ type epsilonWallFunction; Cmu 0.09; kappa 0.41; E 9.8; value $internalField; }}
}}
""")
    scrivi(f"{c}/0/nut", INTESTA.format(cls="volScalarField", obj="nut") + f"""
dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;
boundaryField {{
  inlet {{ type calculated; value uniform 0; }} outlet {{ type calculated; value uniform 0; }}
  lati {{ type calculated; value uniform 0; }} top {{ type calculated; value uniform 0; }}
  ground {{ type nutkAtmRoughWallFunction; Cmu 0.09; kappa 0.41; E 9.8; z0 {z0l}; value uniform 0; }}
}}
""")
    scrivi(f"{c}/0/alphat", INTESTA.format(cls="volScalarField", obj="alphat") + """
dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;
boundaryField {
  inlet { type calculated; value uniform 0; } outlet { type calculated; value uniform 0; }
  lati { type calculated; value uniform 0; } top { type calculated; value uniform 0; }
  ground { type alphatJayatillekeWallFunction; Prt 0.85; value uniform 0; }
}
""")
    scrivi(f"{c}/0/T", INTESTA.format(cls="volScalarField", obj="T") + f"""
dimensions [0 0 0 1 0 0 0]; internalField uniform {TREF};
boundaryField {{
  inlet {{ {prof(th)} }}
  outlet {{ type zeroGradient; }}
  lati {{ type zeroGradient; }} top {{ type fixedValue; value uniform {TREF + Gm * a.alto:.6g}; }}
  ground {{ type zeroGradient; }}
}}
""")
    scrivi(f"{c}/0/p_rgh", INTESTA.format(cls="volScalarField", obj="p_rgh") + f"""
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;
boundaryField {{
  inlet {{ type fixedFluxPressure; value $internalField; }}
  outlet {{ {prof(prgh)} }}
  lati {{ type fixedFluxPressure; value $internalField; }} top {{ type fixedFluxPressure; value $internalField; }}
  ground {{ type fixedFluxPressure; value $internalField; }}
}}
""")
    scrivi(f"{c}/0/p", INTESTA.format(cls="volScalarField", obj="p") + """
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;
boundaryField { inlet { type calculated; value $internalField; } outlet { type calculated; value $internalField; }
  lati { type calculated; value $internalField; } top { type calculated; value $internalField; } ground { type calculated; value $internalField; } }
""")


def stratifica(a):
    """Dopo writeCellCentres: temperatura e pressione iniziali in equilibrio idrostatico (aiuta la convergenza)."""
    c = a.caso
    # maglia di blockMesh a un blocco: punti ordinati (z, y, x), celle con x più veloce, poi y, poi z
    import re
    bm = open(f"{c}/system/blockMeshDict").read()
    nx, ny, nz = (int(v) for v in re.search(r"hex \(0 1 2 3 4 5 6 7\) \((\d+) (\d+) (\d+)\)", bm).groups())
    n, corpo = leggi_lista(f"{c}/constant/polyMesh/points")
    P = np.array(corpo.replace("(", " ").replace(")", " ").split(), float).reshape(nz + 1, ny + 1, nx + 1, 3)
    Z = P[..., 2]
    Cz = 0.125 * (Z[:-1, :-1, :-1] + Z[:-1, :-1, 1:] + Z[:-1, 1:, :-1] + Z[:-1, 1:, 1:]
                  + Z[1:, :-1, :-1] + Z[1:, :-1, 1:] + Z[1:, 1:, :-1] + Z[1:, 1:, 1:]).ravel()
    Gm = gamma(a.N)
    for campo, val in (("T", TREF + Gm * Cz), ("p_rgh", -G0 * BETA * Gm * Cz ** 2 / 2.0)):
        p = f"{c}/0/{campo}"; t = open(p).read()
        vecchio = t[t.index("internalField"):t.index(";", t.index("internalField")) + 1]
        nuovo = f"internalField nonuniform List<scalar> {len(val)}(" + " ".join(f"{v:.7g}" for v in val) + ");"
        open(p, "w").write(t.replace(vecchio, nuovo))
    print(f"campi iniziali stratificati: {len(Cz)} celle, theta da {TREF + Gm * Cz.min():.2f} a {TREF + Gm * Cz.max():.2f} K")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlante"); ap.add_argument("caso")
    ap.add_argument("--sito", default="stretto")
    ap.add_argument("--centro", required=True); ap.add_argument("--da", type=float, required=True)
    ap.add_argument("--N", type=float, required=True, help="frequenza di Brunt-Vaisala, 1/s (0 = neutra)")
    ap.add_argument("--vento", type=float, default=8.0); ap.add_argument("--z0-ingresso", type=float, default=0.0002)
    ap.add_argument("--x", default="-22,22"); ap.add_argument("--y", default="-22,22")
    ap.add_argument("--dx", type=float, default=400.0); ap.add_argument("--alto", type=float, default=10000.0)
    ap.add_argument("--strati", type=int, default=40); ap.add_argument("--primo", type=float, default=10.0)
    ap.add_argument("--sfuma", type=float, default=6.0); ap.add_argument("--liscia", type=float, default=-1.0, help="sigma della lisciatura del rilievo, m (-1 = la maglia)"); ap.add_argument("--delta", type=float, default=500.0); ap.add_argument("--spugna", type=float, default=6000.0)
    ap.add_argument("--iter", type=int, default=2000); ap.add_argument("--processori", type=int, default=4)
    ap.add_argument("--solo-punti", action="store_true"); ap.add_argument("--stratifica", action="store_true")
    a = ap.parse_args()
    a.x = [float(v) for v in a.x.split(",")]; a.y = [float(v) for v in a.y.split(",")]
    if a.liscia < 0:
        a.liscia = a.dx
    lat0, lon0 = (float(v) for v in a.centro.split(","))
    G = Geometria(lat0, lon0, a.da)
    if a.stratifica:
        stratifica(a); return
    A = leggi_atlante(a.atlante)
    if not a.solo_punti:
        nx, ny = caso(a, G)
        print(f"caso scritto: {nx} x {ny} x {a.strati} celle; N {a.N} 1/s, Gamma {gamma(a.N) * 1000:.2f} K/km; ora blockMesh, poi --solo-punti")
        return
    z0 = alza_punti(a, A, G)
    campi(a, z0)
    print(f"punti sul rilievo, {len(z0)} facce di suolo (z0 mediana {np.median(z0):.4f})")


if __name__ == "__main__":
    main()
