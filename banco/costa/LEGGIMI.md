# Vento sotto costa: radar dei Sentinel-1 contro i modelli e Lens (24/11/2026)

Risultati e decisioni: nel Project, `routing/atlante_costa_risultati.md`. Qui gli strumenti e le figure.

- Dati: `shift-grib`, workflow `costa.yml` (`tools/costa.py`), release `costa-sar` (un .npz per area e dataset).
  ECMWF 0,25° dall'archivio di shift-grib (release `archivio-AAAA-MM`).
- `costa_diagnosi.py`: radar ÷ ECMWF e Lens ÷ ECMWF per distanza, ora, stagione, forza (`diagnosi.txt`).
- `costa_ecmwf.py` + `ShiftSar.kt`: ECMWF letto come Shift (normale e solo mare) e Lens sopra, sulle celle del radar.
- `costa_maschere.sh`, `costa_maschere_sintesi.py`: varianti della lettura solo mare.
- `costa_taratura.py`, `costa_pool.py`, `costa_norm.py`: i pezzi di Lens accesi e spenti (proprietà `shift.lens.*`).
- `costa_correzione.py`: prova incrociata di una correzione locale per settore di direzione.
- `costa_figura.py`, `costa_scene.py`, `costa_grafico.py`: le figure (`atlante_<area>.png`, `scene_*.png`,
  `errore_per_distanza.png`).
- `LensSar.kt`: Lens sul vento ECMWF del prodotto (per la diagnosi).

Dati: E.U. Copernicus Marine Service Information (WIND_MED_PHY_HR_L3_NRT_012_104).
