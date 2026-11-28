"""Scarica le immagini Sentinel-1 IW GRDH (VV) dei passaggi di un'area dal Copernicus Data Space, le ritaglia sull'area,
le riduce a celle di 40 m e scrive un piccolo .npz per immagine. Gira su GitHub (segreti CDSE_USER, CDSE_PASS).

  python3 scarica_grd.py <passaggi_area.json> <cartella uscita> <quanti> <gruppo> <gruppi>

Per ogni passaggio (tempo unix): cerca nel catalogo i prodotti IW GRDH che toccano il riquadro entro 10 minuti, scarica
il solo file VV e la sua annotazione (i "nodi" del prodotto), e tiene:
  amp     ampiezza a 40 m, uint16 (radice della media di DN^2 su 4x4), 0 = fuori immagine
  lat,lon latitudine e longitudine ogni 25 celle (1 km), per riportare celle e direzioni sulla carta
  meta    nome del prodotto, tempo, verso dell'orbita
"""
import io, json, os, sys, time, re, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import requests
import tifffile
from scipy.interpolate import RegularGridInterpolator

CAT = "https://catalogue.dataspace.copernicus.eu/odata/v1"
DL = "https://download.dataspace.copernicus.eu/odata/v1"
_tok = {"t": 0, "v": None}


def gettone():
    if time.time() - _tok["t"] > 400:
        r = requests.post("https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
                          data=dict(grant_type="password", client_id="cdse-public", username=os.environ["CDSE_USER"], password=os.environ["CDSE_PASS"]), timeout=60)
        r.raise_for_status(); _tok["v"] = r.json()["access_token"]; _tok["t"] = time.time()
    return {"Authorization": "Bearer " + _tok["v"]}


def prendi(url, stream=False):
    for k in range(4):
        r = requests.get(url, headers=gettone(), timeout=600, stream=stream, allow_redirects=False)
        while r.status_code in (301, 302, 303, 307, 308):
            r = requests.get(r.headers["Location"], headers=gettone(), timeout=600, stream=stream, allow_redirects=False)
        if r.status_code == 200:
            return r
        print("  ", r.status_code, url[-90:], r.text[:200] if not stream else "", flush=True)
        time.sleep(20 * (k + 1))
    return None


def cerca(box, t):
    S, N, W, E = box
    iso = lambda x: time.strftime("%Y-%m-%dT%H:%M:%S.000Z", time.gmtime(x))
    f = ("Collection/Name eq 'SENTINEL-1' and contains(Name,'IW_GRDH') "
         f"and ContentDate/Start lt {iso(t + 600)} and ContentDate/End gt {iso(t - 600)} "
         f"and OData.CSC.Intersects(area=geography'SRID=4326;POLYGON(({W} {S},{E} {S},{E} {N},{W} {N},{W} {S}))')")
    r = requests.get(CAT + "/Products", params={"$filter": f, "$top": 20}, timeout=120)
    r.raise_for_status()
    v = r.json().get("value", [])
    # una sola versione per istante di inizio: meglio quella non COG
    per = {}
    for p in sorted(v, key=lambda p: ("_COG" in p["Name"], p["Name"])):
        per.setdefault(p["ContentDate"]["Start"][:19], p)
    return list(per.values())


def nodi(pid, *percorso):
    url = f"{DL}/Products({pid})" + "".join(f"/Nodes({x})" for x in percorso) + "/Nodes"
    r = prendi(url)
    return [n["Name"] for n in r.json().get("result", [])] if r is not None else []


def file_vv(p, tmp):
    """Scarica il tiff VV e l'annotazione VV. Prima per nodi; se non va, tutto lo zip."""
    pid, nome = p["Id"], p["Name"]
    try:
        mis = [n for n in nodi(pid, nome, "measurement") if "-vv-" in n and n.endswith((".tiff", ".tif"))]
        ann = [n for n in nodi(pid, nome, "annotation") if "-vv-" in n and n.endswith(".xml")]
        if mis and ann:
            x = prendi(f"{DL}/Products({pid})/Nodes({nome})/Nodes(annotation)/Nodes({ann[0]})/$value")
            r = prendi(f"{DL}/Products({pid})/Nodes({nome})/Nodes(measurement)/Nodes({mis[0]})/$value", stream=True)
            if x is not None and r is not None:
                pt = os.path.join(tmp, "vv.tiff")
                with open(pt, "wb") as fh:
                    for ch in r.iter_content(1 << 22):
                        fh.write(ch)
                return pt, x.content
    except Exception as e:  # noqa: BLE001
        print("   nodi:", repr(e)[:200], flush=True)
    print("   scarico tutto lo zip", flush=True)
    r = prendi(f"{DL}/Products({pid})/$value", stream=True)
    if r is None:
        return None, None
    pz = os.path.join(tmp, "p.zip")
    with open(pz, "wb") as fh:
        for ch in r.iter_content(1 << 22):
            fh.write(ch)
    with zipfile.ZipFile(pz) as z:
        mis = [n for n in z.namelist() if "/measurement/" in n and "-vv-" in n]
        ann = [n for n in z.namelist() if re.search(r"/annotation/[^/]*-vv-[^/]*\.xml$", n)]
        pt = os.path.join(tmp, "vv.tiff")
        with open(pt, "wb") as fh:
            fh.write(z.read(mis[0]))
        x = z.read(ann[0])
    os.remove(pz)
    return pt, x


def riduci(p, box, tmp, usc):
    pt, xml = file_vv(p, tmp)
    if pt is None:
        return "scarico fallito"
    root = ET.fromstring(xml)
    pts = [(int(g.find("line").text), int(g.find("pixel").text), float(g.find("latitude").text), float(g.find("longitude").text))
           for g in root.iter("geolocationGridPoint")]
    righe = sorted({q[0] for q in pts}); col = sorted({q[1] for q in pts})
    LA = np.full((len(righe), len(col)), np.nan); LO = LA.copy()
    for li, px, la, lo in pts:
        LA[righe.index(li), col.index(px)] = la; LO[righe.index(li), col.index(px)] = lo
    fla = RegularGridInterpolator((righe, col), LA, bounds_error=False, fill_value=None)
    flo = RegularGridInterpolator((righe, col), LO, bounds_error=False, fill_value=None)
    img = tifffile.imread(pt)
    os.remove(pt)
    ny, nx = img.shape
    rr, cc = np.meshgrid(np.arange(0, ny, 50), np.arange(0, nx, 50), indexing="ij")
    la = fla((rr, cc)); lo = flo((rr, cc))
    S, N, W, E = box
    dentro = (la >= S) & (la <= N) & (lo >= W) & (lo <= E)
    if dentro.sum() < 50:
        return "fuori dal riquadro"
    r0, r1 = rr[dentro].min(), min(ny, rr[dentro].max() + 50); c0, c1 = cc[dentro].min(), min(nx, cc[dentro].max() + 50)
    r1 = r0 + ((r1 - r0) // 4) * 4; c1 = c0 + ((c1 - c0) // 4) * 4
    a = img[r0:r1, c0:c1].astype(np.float32); del img
    vuoto = (a == 0).reshape((r1 - r0) // 4, 4, (c1 - c0) // 4, 4).mean((1, 3)) > 0.25
    amp = np.sqrt((a * a).reshape((r1 - r0) // 4, 4, (c1 - c0) // 4, 4).mean((1, 3)))
    amp = np.clip(np.rint(amp), 1, 65535).astype(np.uint16); amp[vuoto] = 0
    # latitudine e longitudine ogni 25 celle da 40 m (centri delle celle)
    gr, gc = np.meshgrid(np.arange(0, amp.shape[0], 25), np.arange(0, amp.shape[1], 25), indexing="ij")
    lat = fla((r0 + gr * 4 + 1.5, c0 + gc * 4 + 1.5)).astype(np.float32); lon = flo((r0 + gr * 4 + 1.5, c0 + gc * 4 + 1.5)).astype(np.float32)
    orb = root.find(".//pass")
    meta = dict(nome=p["Name"], inizio=p["ContentDate"]["Start"], fine=p["ContentDate"]["End"], orbita=orb.text if orb is not None else "", r0=int(r0), c0=int(c0),
                righe=int(ny), colonne=int(nx), buone=float((amp > 0).mean()))
    np.savez_compressed(usc, amp=amp, lat=lat, lon=lon, meta=json.dumps(meta))
    return f"ok {amp.shape} buone {meta['buone']:.2f} {os.path.getsize(usc) / 1e6:.1f} MB"


if __name__ == "__main__":
    P = json.load(open(sys.argv[1])); usc = sys.argv[2]; quanti = int(sys.argv[3]); gruppo = int(sys.argv[4]); gruppi = int(sys.argv[5])
    os.makedirs(usc, exist_ok=True); tmp = os.environ.get("RUNNER_TEMP", "/tmp")
    fatti = set(os.environ.get("FATTI", "").split())
    for k, t in enumerate(P["t"][:quanti]):
        if k % gruppi != gruppo:
            continue
        try:
            prod = cerca(P["box"], t)
        except Exception as e:  # noqa: BLE001
            print(k, t, "catalogo:", repr(e)[:200], flush=True); continue
        print(k, t, time.strftime("%Y-%m-%d %H:%M", time.gmtime(t)), [p["Name"][:40] for p in prod], flush=True)
        for p in prod:
            nome = f"{P['area']}_{t}_{p['ContentDate']['Start'][11:19].replace(':', '')}.npz"
            if nome in fatti:
                print("   già fatto", nome); continue
            try:
                print("  ", nome, riduci(p, P["box"], tmp, os.path.join(usc, nome)), flush=True)
            except Exception as e:  # noqa: BLE001
                print("   errore:", repr(e)[:300], flush=True)
