# Lens 2: gli analoghi del radar, ECMWF e ICON-2I (25/11/2026)

Segue `lens/studio_statistico_25_11.md`. Giano ha scelto di imparare dal radar area per area. Nell'app il radar non serve: serve per costruire le correzioni. Il metro è il vento vero a 10 m: il radar convertito con la stabilità di ECMWF.

## In breve

- **Il metodo degli analoghi** batte la carta fissa di ieri.
  - Per ogni ora di previsione si cercano, nell'archivio della stessa area, le 25 situazioni passate più simili: vento al largo (forza e direzione), vento a 925 hPa, mattina o sera.
  - Si applica la correzione che il radar aveva misurato in quelle situazioni.
  - Sulle settimane mai viste l'errore sistematico entro 3 km scende del 60-65%.
- **Le correzioni valgono per un modello solo.** Quella imparata su ECMWF, messa sopra ICON-2I, lo peggiora: l'errore sistematico passa da 0,76 a 1,41 nodi.
- **ICON-2I da solo, sulle celle del radar, ha un errore quadratico più alto di ECMWF letto dal mare (4,23 contro 3,89 nodi), e sotto costa sbaglia di più.** Davanti a Palermo, entro 1,5 km dalla costa, sottostima di 4,8 nodi (ECMWF di 3,2), soprattutto con ponente e libeccio. Corretto con i suoi analoghi migliora molto.
- **Il risultato migliore è l'insieme:** ECMWF corretto più ICON-2I corretto, oggi con i pesi 0,6 e 0,4 (vedi il punto 7). Sulle settimane mai viste, con 0,7 e 0,3:
  - errore quadratico −7% rispetto a ECMWF e −15% rispetto a ICON-2I;
  - errore sistematico entro 3 km da 1,43 (ECMWF) e 1,50 (ICON-2I) a 0,47 nodi.
- **youmeteo:** l'archivio ha 10 giorni. Servono 6-12 mesi di passaggi prima di imparare: l'accumulo comincia ora, e resta nel privato.

## Dati

- **Aree:** le 10 di prima più 4 nuove (Pontine, Corsica nord, Monaco, Cicladi), 5.371 passaggi del radar e 6,2 milioni di celle-passaggio.
- **ECMWF:** dati aperti, ora ritagliati su 34,5-45,5 N e 5-27 E. Vento a 10 m e campi di stabilità, giugno 2024 - novembre 2026.
- **ICON-2I:** la storia di Open-Meteo, dall'aprile 2025: 3.576 ore attorno ai passaggi (6 al giorno) per 13 aree (le Cicladi sono fuori dal suo dominio).
  - Letto come Shift: bilineare dai soli punti di mare della sua griglia (0,02° × 0,025°).
  - 3.203 passaggi in comune con ECMWF: il confronto è fatto sulle stesse celle e sulle stesse ore.
- **Prove:**
  - **settimane:** si impara sulle settimane pari e si giudica sulle dispari, per scegliere i parametri;
  - **cassaforte:** una settimana su cinque, tenuta da parte fino alla prova finale.

## 1. I metodi area per area (ECMWF, 14 aree, cassaforte)

| Metodo | Errore quadratico | Sistematico | Entro 3 km | 3-10 km |
|---|---|---|---|---|
| ECMWF letto dal mare | 3,92 | 1,03 | 1,60 | 1,06 |
| Carta fissa (cella × ottante) | 3,81 | 0,66 | 0,80 | 0,63 |
| Armoniche (per cella, continua nella direzione) | 3,78 | 0,56 | 0,77 | 0,50 |
| Analoghi (vento al largo, stabilità, ora) | 3,75 | 0,44 | 0,61 | 0,40 |
| **Analoghi (vento al largo, 925 hPa, ora)** | **3,75** | **0,42** | **0,56** | **0,34** |

Nodi, vento vero.
- **L'errore quadratico scende poco perché è fatto quasi tutto di rumore del radar** (3-4 nodi su una cella) e di errori di tempo del modello.
- **La variante con il vento a 925 hPa al posto della stabilità** va bene uguale, e nell'app è più comoda. Il vento a 925 hPa c'è già nei file di profilo di shift-grib; per la stabilità servirebbe anche la temperatura del mare.

Entro 3 km, per area (cassaforte, nodi): modello → analoghi.

| Area | Modello | Analoghi |
|---|---|---|
| Palermo | 2,55 | 1,64 |
| Messina | 2,81 | 0,00 (sotto il rumore) |
| Genova | 2,06 | 0,00 |
| Cagliari | 2,70 | 1,17 |
| Egadi | 1,25 | 0,22 |
| Eolie | 2,08 | 1,50 |
| Bonifacio e La Maddalena | 1,02 | 0,31 |
| Pontine | 1,47 | 0,41 |
| Corsica nord | 1,62 | 0,00 |
| Cicladi | 1,56 | 0,48 |

"0,00" vuol dire che l'errore sistematico che resta è più piccolo di quanto il rumore del radar lasci misurare.

**Quanti passaggi servono** (prova sulle settimane, imparando su una frazione dei passaggi):

| Passaggi per imparare, per area | Errore sistematico | Entro 3 km |
|---|---|---|
| nessuno (modello) | 1,07 | 1,73 |
| circa 15 | 1,26 | 1,53 |
| circa 35 | 1,04 | 1,33 |
| circa 70 | 0,84 | 1,11 |
| circa 140 | 0,72 | 0,95 |

- **Sotto i 40 passaggi non conviene.** Da 70 in su sì.
- **Con più passaggi migliora ancora:** nella prova finale si impara su circa l'80% (200-400 per area).
- Il radar passa su un'area 10-15 volte al mese, quindi le correzioni migliorano da sole aggiornandole ogni mese.

**Provati e non utili** (dentro il rumore):
- la forma del vento del modello sull'area (componenti principali);
- la stagione;
- il contrasto terra-mare;
- più o meno analoghi tra 15 e 60.

Lens 0.54 sotto gli analoghi dà un piccolo guadagno: va tenuto come base, perché disegna anche sotto il chilometro, dove il radar non arriva.

## 2. ECMWF e ICON-2I (13 aree, stesse celle e stesse ore)

Cassaforte:

| | Errore quadratico | Sistematico | Entro 3 km | Errore medio |
|---|---|---|---|---|
| ECMWF letto dal mare | 3,89 | 0,95 | 1,43 | −0,15 |
| ICON-2I | 4,23 | 0,76 | 1,50 | +0,25 |
| Media dei due, senza correzione (0,65/0,35) | 3,67 | 0,71 | 1,21 | −0,01 |
| ECMWF + analoghi | 3,78 | 0,61 | 0,59 | +0,06 |
| ICON-2I + analoghi | 4,07 | 0,63 | 0,86 | −0,34 |
| ICON-2I con la correzione di ECMWF | 4,67 | 1,41 | 1,65 | +0,57 |
| **Insieme corretto: 0,7 ECMWF + 0,3 ICON-2I** | **3,61** | **0,52** | **0,47** | −0,06 |

Sulle settimane la classifica è la stessa (insieme corretto 2,94 contro 3,30 di ECMWF e 3,44 di ICON-2I). Il peso migliore sta tra 0,6 e 0,7 per ECMWF.

**Perché ICON-2I da solo va peggio:**
- **I suoi dettagli spesso sono fuori posto o fuori tempo.** Un getto spostato di qualche chilometro conta due volte nell'errore: dove c'è e non dovrebbe, e dove manca. ECMWF, più liscio, sbaglia meno così.
- **Vicino alla costa sente troppo la terra.** Palermo con ponente e libeccio, entro 3 km: ICON-2I −3,3 nodi, ECMWF −0,4/−0,7. Nelle scene della cassaforte il golfo di ICON-2I è spesso in calma dove il radar vede 20-25 nodi (figura `banco/figure/modelli_cassaforte.png`).
- **Il radar passa alle 05 e alle 17 UTC.** Le brezze del pomeriggio, dove ICON-2I dovrebbe rendere di più, non sono misurate.

**Risposta alla domanda di Giano** ("tarando su ECMWF, gli altri modelli migliorano, peggiorano o restano uguali?"):
- **Peggiorano,** se gli si applica la correzione di ECMWF.
- **Migliorano,** se ognuno ha la sua.
- **Insieme danno il risultato migliore.**

## 3. Lens 2, come entrerebbe nell'app

- **Per area e per modello, una "libreria":**
  - per ogni passaggio imparato, i dati della situazione: vento al largo, vento a 925 hPa, mattina o sera;
  - la carta del radar ÷ modello a 1 km, più cinque carte delle armoniche come base;
  - due byte per cella e passaggio: Palermo circa 1 MB, l'Arcipelago toscano 10-20 MB. Si scarica per area, come gli atlanti di Lens.
- **A ogni ora di previsione:**
  1. si calcolano i dati della situazione dal GRIB;
  2. si cercano i 25 passaggi più simili;
  3. si fa la media pesata delle loro carte, tirata verso le armoniche dove i dati sono pochi;
  4. il fattore si applica sopra Lens 0.54 (che resta per i dettagli sotto il chilometro).
- **Con più modelli caricati,** l'insieme corretto: 0,6 ECMWF + 0,4 ICON-2I (punto 7). Poi, quando avrà abbastanza storia, youmeteo.
- **La direzione non si corregge.** Il prodotto del radar prende la direzione da ECMWF, quindi non la misura davvero.
- **Aggiornamento:** una volta al mese si aggiungono i passaggi nuovi. È lavoro del deposito `lens` e non costa niente.

## 4. Bonacce: un primo sguardo

Estate, vento al largo sotto 4 m/s, celle entro 25 km dalla costa. La calma qui è sotto 4 nodi.
- **Il modello mette più calme del radar:** a Palermo il 49% delle celle contro il 31%, alle Eolie il 56% contro il 39%.
- **Dove il modello dà calma,** il radar vede vento nel 20-25% delle celle. Il contrario succede nell'8-10%.

**Con due cautele forti:**
- sotto 2 m/s il radar è poco affidabile;
- d'estate, con aria stabile, il vento vero è più forte di quello neutro che il radar misura.

I grandi buchi d'aria del largo, quelli dell'alta stazionaria, vanno studiati a parte. Servono riquadri al largo lungo le rotte (per esempio Palermo-Bonifacio-Monaco) e ICON-2I su tutte le celle, non solo dove c'è vento. È il prossimo studio, se interessa.

## 5. Limiti

- **Radar:** passa alle 05 e alle 17 UTC, ha 3-4 nodi di rumore per cella, e nel primo chilometro dalla costa può avere difetti suoi.
- **ICON-2I** è quello della storia di Open-Meteo: per ogni ora le prime ore della corsa più recente. Con anticipi di 1-3 giorni gli errori del modello crescono. La correzione della costa dovrebbe reggere, ma va verificato con l'archivio di shift-grib quando sarà lungo abbastanza.
- **Le stime per area nella cassaforte sono rumorose:** in qualche area ci sono pochi gruppi di passaggi.
- **Monaco:** troppe poche celle entro 3 km con abbastanza passaggi per la misura vicino alla costa.

## 6. Prossimi passi

1. **Il motore in Kotlin** in `lens-motore`: la libreria degli analoghi, la ricerca, l'insieme dei modelli. Casi di prova presi da questi numeri.
2. **Le librerie** per 14 aree con ECMWF e 13 con ICON-2I, nel formato per l'app.
3. **youmeteo:** l'accumulo dall'archivio del deposito privato. La prima prova quando ci saranno circa 70 passaggi per area (primavera).
4. **Bonacce lungo le rotte:** riquadri al largo e ICON-2I completo.
5. **Altre aree a richiesta:** Corsica ovest, Liguria di ponente, Sardegna nord-est.

## 7. Verifiche prima della consegna (25/11, pomeriggio)

Fatte su richiesta di Giano: prima di passare Lens 2 a Shift, che cosa si può ancora migliorare e che cosa va controllato.

**Previsioni a uno e due giorni.** Le librerie sono imparate su corse con 5 ore di anticipo. Ho scaricato le corse ECMWF con 29 e 53 ore di anticipo per tutti i 5.371 passaggi, con il loro vento a 925 hPa. Cassaforte, 14 aree, nodi:

| Anticipo | Libreria imparata a | Errore quadratico (modello → corretto) | Sistematico | Entro 3 km |
|---|---|---|---|---|
| 5 ore | 5 ore | 3,92 → 3,74 | 1,03 → 0,42 | 1,60 → 0,57 |
| 29 ore | 5 ore | 4,08 → 3,94 | 1,10 → 0,51 | 1,57 → 0,69 |
| 29 ore | 29 ore | 4,08 → 3,94 | 1,10 → 0,56 | 1,57 → 0,70 |
| 53 ore | 5 ore | 4,28 → 4,21 | 1,06 → 0,50 | 1,48 → 0,58 |
| 53 ore | 53 ore | 4,28 → 4,20 | 1,06 → 0,50 | 1,48 → 0,53 |

- **La correzione regge a uno e due giorni.** L'errore sistematico vicino alla costa si dimezza o meno ancora.
- **L'errore quadratico scende meno** (−3% a un giorno, −2% a due): a quegli anticipi pesano gli errori di tempo del modello, che nessuna correzione della costa può togliere.
- **Librerie per anticipo non servono:** quella imparata a 5 ore va bene uguale. Si tiene una libreria sola per area e modello.

**Le boe, a tutte le ore.** Nel prodotto in situ di Copernicus Marine ci sono solo due boe con il vento dentro o vicino alle nostre aree:
- **Ponza** (marzo-ottobre 2026, a 1,7 km dall'isola);
- **Calvi** (Météo-France, marzo 2025 - settembre 2026, al largo: la cella più vicina con il radar è a 10 km).

Per allargare la prova ho esteso la Corsica nord fino a 8,25° E. Confronti ogni 3 ore (corse delle 00 e 12, passi 3-12 ore). Il vento della boa è a pochi metri sul mare, quindi ogni previsione è moltiplicata per un solo numero stimato su tutta la serie. Errore quadratico in nodi:

| | Ponza, tutte | Ponza, 12-15 UTC | Calvi, tutte | Calvi, 12-15 UTC |
|---|---|---|---|---|
| ECMWF | 2,76 | 2,34 | 3,24 | 3,07 |
| ECMWF corretto | 2,80 | 2,39 | 3,21 | 3,03 |
| ICON-2I | 2,98 | 2,41 | 3,51 | 3,15 |
| ICON-2I corretto | 3,00 | 2,41 | 3,52 | 3,18 |
| Insieme corretto (0,6/0,4) | **2,51** | **2,06** | **2,98** | **2,73** |

- **Nei due punti delle boe la correzione è quasi neutra:** lì il radar vede poco da correggere (fattore 1,01-1,08), e la boa di Ponza sente il ridosso dell'isola più fine del chilometro del radar (con la tramontana 0,88 contro 1,01).
- **Non peggiora a mezzogiorno,** le ore che il radar non vede.
- **L'insieme dei due modelli migliora molto anche qui, a tutte le ore:** −9% a Ponza, −8% a Calvi. È il guadagno più robusto.
- **Le due boe non bastano a provare le correzioni grandi** (golfo di Palermo, Stretto, Genova): lì le boe non ci sono. Le stazioni a terra sui moli sono il passo successivo.

**La fiducia.** La dispersione degli analoghi attorno al fattore dice quanto fidarsi. Sulle settimane di prova, divise in cinque parti uguali:

| Dispersione | Errore relativo che resta |
|---|---|
| sotto 0,17 | 0,18 |
| 0,17-0,22 | 0,23 |
| 0,22-0,26 | 0,26 |
| 0,26-0,33 | 0,30 |
| sopra 0,33 | 0,38 |

È nel motore (`fattoreEDispersione`): Shift la può mostrare come carta della fiducia o usare nel routing.

**Provati e non utili:**
- **analoghi per zone di 15-40 km nelle aree grandi:** uguali agli analoghi sull'area intera;
- **pesi dell'insieme diversi per area o per distanza dalla costa:** il peso migliore di ECMWF sta tra 0,55 e 0,65 dappertutto.

Si usa 0,6 e 0,4 invece di 0,7 e 0,3: l'errore quadratico è un po' più basso (cassaforte 3,60 contro 3,61; settimane 2,92 contro 2,94).

**Un errore da non fare, visto in prova:** una libreria applicata al modello sbagliato rovina la correzione. Nel motore c'è `valePer(modello)`.

## File

- **Codice in `banco/strumenti/`:**
  - `metodi_area.py`: carta, armoniche, analoghi;
  - `banco_icon2i.py`: ICON-2I sulle celle del radar;
  - `figura_scene.py`, `figura_modelli.py`: le figure.
- **Dati in `banco/dati/`:**
  - `icon2i_storia.py`: la storia di ICON-2I da Open-Meteo, flusso `icon2i.yml`;
  - `ecmwf_aperto.py`: ritaglio a scelta;
  - aree nuove in `costa.py`.
- **Figure in `banco/figure/`:**
  - `scene_cassaforte.png`: radar, ECMWF, ECMWF con gli analoghi;
  - `modelli_cassaforte.png`: radar, ECMWF, ICON-2I, insieme corretto.
- **Numeri in `banco/risultati/`:**
  - `analoghi_14_aree_*.json`;
  - `modelli_13_aree.json`.
