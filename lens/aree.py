"""Le aree pronte di Lens: le stesse di app/src/main/java/app/regata/core/Lens.kt (Lens.AREAS). Sud, nord, ovest, est.
Se si cambia un'area qui, va cambiata anche li' (e l'atlante va ricostruito)."""
AREE = {
    # Giano 23/11 sera: Egadi con lo Stagnone, Eolie fino alla costa, Malta piu' larga, Toscano piu' su, Minorca fino al
    # canale, Cicladi piu' larghe, dal Pireo a Kea, Canarie
    "egadi": (37.70, 38.15, 11.95, 12.60),
    "eolie": (38.00, 38.86, 14.25, 15.30),
    "palermo": (37.945, 38.287, 13.108, 13.612),
    "messina": (37.90, 38.35, 15.40, 15.75),
    "malta": (35.68, 36.20, 14.05, 14.75),
    "napoli": (40.50, 40.90, 13.80, 14.62),
    "cagliari": (38.85, 39.25, 8.85, 9.65),
    "maddalena-bonifacio": (41.00, 41.50, 8.95, 9.75),
    "toscano": (42.20, 43.35, 9.75, 11.30),
    "genova": (44.00, 44.47, 8.40, 9.50),
    # Giano, 25/11: le aree di Lens 2 per una Palermo-Monte Carlo, con gli stessi riquadri delle librerie del radar
    "pontine": (40.70, 41.32, 12.70, 13.65),
    "corsica-nord": (42.45, 43.10, 8.25, 9.65),
    "monaco": (43.45, 43.85, 7.15, 8.00),
    # Giano, 21/11
    "venezia": (45.15, 45.60, 12.10, 12.75),
    "kvarner": (44.40, 45.35, 14.00, 15.05),
    "kornati": (43.55, 44.30, 14.85, 15.75),
    "dalmazia": (42.85, 43.55, 15.95, 17.25),
    "maiorca": (39.05, 40.03, 2.22, 3.57),
    "minorca": (39.70, 40.15, 3.50, 4.45),
    "ibiza": (38.55, 39.18, 1.08, 1.78),
    "pireo": (37.45, 38.00, 23.55, 24.55),
    "cicladi": (36.70, 37.75, 24.40, 25.80),
    "saronico": (37.15, 37.95, 22.75, 23.75),
    "capi": (36.10, 36.90, 22.20, 23.30),
    "ionio": (37.95, 38.95, 20.30, 21.20),
    "canarie": (27.70, 28.65, -16.95, -15.30),
}

# 0.57 (Giano, 25/11): le aree fini, a 100 m, dentro le aree grandi (Lens 2 usa le librerie dell'area grande: PADRE).
FINI = {
    "egadi-favignana": (37.85, 38.10, 12.20, 12.52),
    "egadi-marettimo": (37.88, 38.06, 11.95, 12.17),
    "eolie-lipari": (38.35, 38.60, 14.80, 15.04),
    "eolie-stromboli": (38.58, 38.86, 15.00, 15.30),
    "eolie-filicudi": (38.45, 38.68, 14.28, 14.65),
    "palermo-fine": (38.04, 38.27, 13.15, 13.58),
    "napoli-capri": (40.50, 40.70, 14.10, 14.45),
    "napoli-ischia": (40.65, 40.85, 13.82, 14.15),
    "napoli-golfo": (40.70, 40.88, 14.12, 14.45),
    "bonifacio-bocche": (41.20, 41.45, 9.00, 9.32),
    "bonifacio-maddalena": (41.08, 41.33, 9.25, 9.60),
    "pontine-ponza": (40.82, 41.02, 12.80, 13.12),
    "pontine-ventotene": (40.72, 40.88, 13.33, 13.55),
    "pontine-circeo": (41.13, 41.32, 12.95, 13.25),
    "pontine-gaeta": (41.10, 41.30, 13.40, 13.68),
}
PADRE = {"egadi-favignana": "egadi", "egadi-marettimo": "egadi", "eolie-lipari": "eolie", "eolie-stromboli": "eolie", "eolie-filicudi": "eolie", "palermo-fine": "palermo", "napoli-capri": "napoli", "napoli-ischia": "napoli", "napoli-golfo": "napoli", "bonifacio-bocche": "maddalena-bonifacio", "bonifacio-maddalena": "maddalena-bonifacio", "pontine-ponza": "pontine", "pontine-ventotene": "pontine", "pontine-circeo": "pontine", "pontine-gaeta": "pontine"}
MAGLIA_FINI = 100
AREE.update(FINI)

if __name__ == "__main__":
    # per il workflow: l'elenco delle aree chieste come matrice JSON ("tutte", oppure id separati da virgole)
    import json, os, sys
    chieste = (os.environ.get("AREE") or "tutte").replace(" ", "")
    if chieste == "fini":
        chieste = ",".join(FINI)
    box = (os.environ.get("AREA") or "").strip()
    if box:
        ident = "".join(c for c in (os.environ.get("ID") or "area") if c.isalnum() or c == "-")[:40] or "area"
        voci = [{"id": ident, "box": box, "dx": ""}]
    else:
        ids = list(AREE) if chieste in ("", "tutte") else [x for x in chieste.split(",") if x in AREE]
        voci = [{"id": k, "box": ",".join(str(x) for x in AREE[k]), "dx": str(MAGLIA_FINI) if k in FINI else ""} for k in ids]
    if not voci:
        sys.exit("nessuna area")
    print(json.dumps({"include": voci}))
