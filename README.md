# Lens

Laboratorio per il vento sotto costa ad alta definizione: come il rilievo, la costa e il caldo del giorno cambiano il vento di un modello meteo. Il risultato finale è un motore leggero che gira dentro un'app di navigazione, senza rete, sul GRIB già scaricato.

## Che cosa c'è

- `lens/`: il primo Lens (novembre 2026). Ridistribuisce il vento del modello sul rilievo vero, a 150-300 m, con la conservazione della massa, e aggiunge ridossi e raffiche. Il metodo è in `lens/METODO.md`.
- `banco/costa/`: il confronto con il vento misurato dal radar dei satelliti Sentinel-1, sulle dieci aree di prova (24/11/2026).
- `banco/strumenti/`: controlli e prove del laboratorio.

## Regole

- Si usano solo dati aperti, e se ne citano le fonti.
- Ogni modifica si giudica sui dati, con la stessa regola per tutti:
  - si tara su una parte dei giorni e si verifica sugli altri;
  - una parte finale non si guarda fino alla fine.
- Resta solo quello che migliora sui giorni di verifica.

## Dati e licenze

- Vento dal radar: E.U. Copernicus Marine Service Information, prodotto `WIND_MED_PHY_HR_L3_NRT_012_104`.
- ECMWF: dati aperti, CC BY 4.0.
- ICON-2I e ICON-2I-RUC: dati aperti ItaliaMeteo-Arpae (MeteoHub).
- Rilievo: Copernicus DEM. Uso del suolo: ESA WorldCover.
