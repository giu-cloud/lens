"""Corrente di superficie ogni ora (Copernicus Marine, Mediterraneo, con le maree) su un riquadro: per lo Stretto (26/11).
  python correnti.py --box 37.9,38.35,15.4,15.75 --da 2024-04-01 --a 2026-11-25 --out correnti
Dati: E.U. Copernicus Marine Service Information, MEDSEA_ANALYSISFORECAST_PHY_006_013 (NEMO 1/24°, con le maree).
"""
import argparse, os, sys, traceback
import numpy as np

DATASET = "cmems_mod_med_phy-cur_anfc_4.2km-2D_PT1H-m"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--box", default="37.9,38.35,15.4,15.75")
    ap.add_argument("--da", default="2024-04-01"); ap.add_argument("--a", default="2026-11-25")
    ap.add_argument("--out", default="correnti"); ap.add_argument("--nome", default="stretto")
    o = ap.parse_args()
    os.makedirs(o.out, exist_ok=True)
    s, n, w, e = [float(x) for x in o.box.split(",")]
    import copernicusmarine as cm
    import xarray as xr
    fn = f"raw_{o.nome}.nc"
    try:
        cm.subset(dataset_id=DATASET, variables=["uo", "vo"], minimum_longitude=w, maximum_longitude=e,
                  minimum_latitude=s, maximum_latitude=n, start_datetime=f"{o.da}T00:00:00", end_datetime=f"{o.a}T23:59:59",
                  output_directory=o.out, output_filename=fn, username=os.environ.get("CMEMS_USER"),
                  password=os.environ.get("CMEMS_PASS"), disable_progress_bar=True)
    except Exception:
        traceback.print_exc()
        try:
            print(cm.describe(dataset_id=DATASET))
        except Exception:
            traceback.print_exc()
        return 1
    ds = xr.open_dataset(os.path.join(o.out, fn))
    print(ds)
    u = ds["uo"].values; v = ds["vo"].values
    while u.ndim > 3:
        u = u[:, 0]; v = v[:, 0]
    t = ds["time"].values.astype("datetime64[s]").astype("int64")
    np.savez_compressed(os.path.join(o.out, f"correnti_{o.nome}.npz"), t=t, lat=ds["latitude"].values, lon=ds["longitude"].values,
                        u=u.astype("f4"), v=v.astype("f4"))
    os.remove(os.path.join(o.out, fn))
    print("fatto", u.shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
