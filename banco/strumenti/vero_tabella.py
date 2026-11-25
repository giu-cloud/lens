"""Il radar della tabella in vento vero a 10 m (vento_vero.neutro_a_vero con la stabilità di ECMWF della scena).

  python3 vero_tabella.py <tabella.npz> <uscita.npy>

Dove la stabilità manca, resta il vento neutro del radar.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vento_vero import neutro_a_vero


def main():
    tab, usc = sys.argv[1:3]
    T = np.load(tab)
    r = T["radar"].astype("f8")
    ok = np.isfinite(r) & np.isfinite(T["dTam"]) & np.isfinite(T["sst"]) & np.isfinite(T["umid"])
    v = r.copy()
    v[ok] = neutro_a_vero(r[ok], T["dTam"][ok], T["sst"][ok], T["umid"][ok])
    np.save(usc, v.astype("f4"))
    print(f"{ok.mean() * 100:.1f}% delle righe con la stabilità; vero ÷ neutro mediano {np.nanmedian(v[ok] / r[ok]):.3f}")


if __name__ == "__main__":
    main()
