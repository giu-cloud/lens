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
    "eolie-ovest": (38.44, 38.67, 14.25, 14.7),
    "eolie-centro": (38.29, 38.7, 14.72, 15.12),
    "eolie-est": (38.55, 38.86, 14.97, 15.3),
    "eolie-milazzo": (38.15, 38.35, 15.05, 15.3),
    "palermo-fine": (38.0, 38.287, 13.108, 13.612),
    "egadi-ovest": (37.86, 38.1, 11.95, 12.22),
    "egadi-est": (37.75, 38.12, 12.18, 12.58),
    "napoli-ovest": (40.62, 40.9, 13.8, 14.25),
    "napoli-est": (40.5, 40.86, 14.15, 14.55),
    "bonifacio-ovest": (41.15, 41.5, 8.95, 9.3),
    "bonifacio-est": (41.05, 41.42, 9.25, 9.7),
    "pontine-ponza": (40.8, 41.05, 12.75, 13.15),
    "pontine-ventotene": (40.7, 40.9, 13.3, 13.58),
    "pontine-circeo": (41.1, 41.32, 12.9, 13.3),
    "pontine-gaeta": (41.1, 41.32, 13.25, 13.65),
    "toscano-capraia": (42.93, 43.17, 9.75, 10.0),
    "toscano-elba": (42.66, 42.98, 10.02, 10.58),
    "toscano-giglio": (42.2, 42.5, 10.78, 11.3),
    "malta-gozo": (35.88, 36.17, 14.08, 14.5),
    "malta-valletta": (35.72, 35.96, 14.25, 14.68),
    "corsica-capocorso": (42.75, 43.1, 9.18, 9.62),
    "maiorca-pollenca": (39.75, 40.03, 2.95, 3.35),
    "maiorca-tramuntana": (39.7, 39.95, 2.5, 2.95),
    "maiorca-dragonera": (39.5, 39.75, 2.22, 2.6),
    "maiorca-palma": (39.35, 39.6, 2.45, 2.85),
    "dalmazia-split": (43.25, 43.55, 16.1, 16.7),
    "dalmazia-hvar": (43.05, 43.32, 16.3, 16.9),
    "dalmazia-vis": (42.95, 43.15, 15.95, 16.3),
    "dalmazia-korcula": (42.85, 43.05, 16.55, 17.25),
    "canarie-tenerife": (27.95, 28.65, -16.95, -16.1),
    "canarie-grancanaria": (27.7, 28.25, -15.92, -15.3),
    "cagliari-golfo": (38.95, 39.25, 9.0, 9.42),
    "cagliari-pula": (38.85, 39.1, 8.85, 9.05),
    "cagliari-carbonara": (38.98, 39.22, 9.4, 9.65),
    "genova-ovest": (44.15, 44.45, 8.4, 8.85),
    "genova-est": (44.2, 44.47, 8.85, 9.5),
    "monaco-nizza": (43.55, 43.8, 7.15, 7.5),
    "monaco-sanremo": (43.65, 43.85, 7.45, 7.9),
    "cicladi-tinos-mykonos": (37.38, 37.66, 25.0, 25.48),
    "cicladi-syros": (37.33, 37.55, 24.8, 25.05),
    "cicladi-paros-naxos": (36.95, 37.2, 25.0, 25.6),
}
# 0.57.1 (Giano, 25/11 sera): allargate, circa 10 km di mare a 100 m attorno a ogni gruppo di isole, dentro l'area madre.
# 25/11 notte (Giano): Arcipelago toscano (Capraia, Elba e canale di Piombino, Giglio e Argentario) e Malta (Gozo e Comino, Malta).
# 26/11 (Giano): Capo Corso, Maiorca, Dalmazia, Canarie, Cagliari, Genova, Nizza-Monaco, Cicladi.
PADRE = {"eolie-ovest": "eolie", "eolie-centro": "eolie", "eolie-est": "eolie", "eolie-milazzo": "eolie", "palermo-fine": "palermo", "egadi-ovest": "egadi", "egadi-est": "egadi", "napoli-ovest": "napoli", "napoli-est": "napoli", "bonifacio-ovest": "maddalena-bonifacio", "bonifacio-est": "maddalena-bonifacio", "pontine-ponza": "pontine", "pontine-ventotene": "pontine", "pontine-circeo": "pontine", "pontine-gaeta": "pontine", "toscano-capraia": "toscano", "toscano-elba": "toscano", "toscano-giglio": "toscano", "malta-gozo": "malta", "malta-valletta": "malta", "corsica-capocorso": "corsica-nord", "maiorca-pollenca": "maiorca", "maiorca-tramuntana": "maiorca", "maiorca-dragonera": "maiorca", "maiorca-palma": "maiorca", "dalmazia-split": "dalmazia", "dalmazia-hvar": "dalmazia", "dalmazia-vis": "dalmazia", "dalmazia-korcula": "dalmazia", "canarie-tenerife": "canarie", "canarie-grancanaria": "canarie", "cagliari-golfo": "cagliari", "cagliari-pula": "cagliari", "cagliari-carbonara": "cagliari", "genova-ovest": "genova", "genova-est": "genova", "monaco-nizza": "monaco", "monaco-sanremo": "monaco", "cicladi-tinos-mykonos": "cicladi", "cicladi-syros": "cicladi", "cicladi-paros-naxos": "cicladi"}
MAGLIA_FINI = 100
# le isole grandi delle Canarie non stanno in 165 mila celle a 100 m
MAGLIA = {"canarie-tenerife": 200, "canarie-grancanaria": 150}
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
        voci = [{"id": k, "box": ",".join(str(x) for x in AREE[k]), "dx": str(MAGLIA.get(k, MAGLIA_FINI)) if k in FINI else ""} for k in ids]
    if not voci:
        sys.exit("nessuna area")
    print(json.dumps({"include": voci}))
