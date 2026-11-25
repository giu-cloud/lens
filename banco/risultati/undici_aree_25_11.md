# Lens 2 su undici aree nuove, e un Lens 2 per un'area qualsiasi (25/11/2026)

Giano, 25/11: "lavora pure a tutte le altre aree a cui riesci a lavorare", e: un Lens ad hoc per un'area scelta a piacere.

## 1. Le aree

Sono le aree di Lens nell'app che non avevano ancora Lens 2, con gli stessi riquadri:
- Adriatico: Venezia, Quarnaro (`kvarner`), Incoronate (`kornati`), Dalmazia;
- Baleari: Maiorca, Minorca, Ibiza;
- Grecia: Pireo, Saronico, capi del Peloponneso (`capi`), Ionio.

Venti e rilievi nuovi rispetto alle 14 aree di prima: bora, maestrale delle Baleari, meltemi dell'Egeo.

- **Radar:** da 260 (Venezia) a 512 (Maiorca) passaggi dal 2024.
- **ECMWF:** scaricato in più per le Baleari e Venezia (riquadro 38,25-45,85 N, 0,75-13,25 E).
- **ICON-2I:** dove arriva la sua griglia, cioè le aree adriatiche, lo Ionio, Minorca e metà di Maiorca. Le aree greche ne restano fuori.
- **Le Canarie restano fuori:** il prodotto radar che uso copre solo il Mediterraneo.

## 2. Quanto rende: la cassaforte (settimane mai usate per imparare)

ECMWF, nodi: errore quadratico, sistematico, sistematico entro 3 km dalla costa. Modello → corretto.

| Area | Errore quadratico | Sistematico | Entro 3 km |
|---|---|---|---|
| Quarnaro | 5,56 → 4,65 | 2,66 → 0,78 | 3,27 → 0,98 |
| Dalmazia | 4,55 → 3,71 | 2,58 → 0,19 | 4,04 → 0,69 |
| Incoronate | 3,55 → 3,50 | 0,77 → 0,74 | 0,77 → 0,15 |
| Ionio | 3,62 → 3,33 | 1,39 → 0,65 | 1,70 → 0,27 |
| Pireo | 3,51 → 3,08 | 1,62 → 0,72 | 1,77 → 0,81 |
| Saronico | 3,49 → 3,09 | 1,41 → 0,63 | 1,63 → 0,63 |
| capi del Peloponneso | 4,06 → 3,69 | 1,27 → 0,32 | 1,29 → 0,45 |
| Venezia | 3,74 → 3,66 | 0,65 → 0,54 | 1,17 → meno del rumore |
| Maiorca | 3,20 → 3,12 | 0,81 → 0,50 | 1,51 → 0,80 |
| Minorca | 3,48 → 3,30 | 1,19 → 0,51 | 2,10 → 0,45 |
| Ibiza | 2,83 → 2,59 | 0,91 → 0,76 | 1,45 → 1,00 |
| **Le sette di Adriatico e Grecia insieme** | **4,01 → 3,54** | **1,66 → 0,61** | **2,21 → 0,65** |
| **Le quattro di Venezia e Baleari insieme** | **3,24 → 3,10** | **0,97 → 0,57** | **1,68 → 0,75** |

- **Il metodo regge lontano da dove l'ho messo a punto.** Nelle aree della bora rende più che altrove. Sul Quarnaro e in Dalmazia l'errore quadratico scende del 16-19% e l'errore sistematico quasi sparisce. È la bora che esce accelerata dai valichi e dai canali tra le isole: ECMWF a 25 km non la vede.
- **Le aree già "facili" migliorano poco.** Le Incoronate al largo e Venezia avevano poco errore sistematico da togliere.

**ICON-2I,** errore quadratico, modello → corretto:

| Area | ICON-2I |
|---|---|
| Quarnaro | 4,35 → 4,16 |
| Incoronate | 4,25 → 4,06 |
| Dalmazia | 3,63 → 3,50 |
| Ionio | 4,01 → 3,61 |
| Minorca | 3,35 → 3,29 |
| Venezia | 2,77 → 2,82 (peggiora) |
| Maiorca | 3,52 → 3,58 (peggiora) |

**Venezia e Maiorca restano senza libreria ICON-2I,** perché lì la correzione peggiora:
- **Venezia:** pochi passaggi del radar (98 per imparare), e ICON-2I già da solo è il migliore. A Venezia conviene usare ICON-2I senza correzione: meglio anche della media dei due modelli (2,77 contro 2,94).
- **Maiorca:** ICON-2I copre solo la metà est dell'isola.

**La media dei due modelli corretti, ECMWF e ICON-2I,** sulle quattro aree adriatiche e ioniche:
- media 0,6/0,4: errore quadratico 3,24;
- ECMWF corretto da solo: 3,50;
- ICON-2I corretto da solo: 3,79.

Il peso migliore sta tra 0,5 e 0,7, come nelle altre aree: resta 0,6/0,4.

## 3. Pubblicato

- **Le librerie:** release `librerie-lra1`, stesso indirizzo. Ora sono 43: 25 aree con ECMWF, 18 anche con ICON-2I.
- **L'indice è rifatto** con il campo `attiva` (tutte true) e una nota su Venezia e Maiorca.
- **L'app le trova da sola** al controllo giornaliero dell'indice. Le aree sono già in `Lens.AREAS` con i loro atlanti, quindi Lens 2 si applica sopra Lens come nelle altre.

## 4. Un Lens 2 per un'area qualsiasi

**Si può fare, ed è quello che ho fatto oggi a mano undici volte.** Per un riquadro qualsiasi nel Mediterraneo:

| Passo | Tempo sulle macchine di GitHub |
|---|---|
| il radar dal 2024 (Copernicus Marine, account già nei segreti) | circa 5 minuti |
| ECMWF alle ore del radar, solo se il riquadro esce da quello già scaricato (ora 34,5-45,85 N, 0,75-27 E) | circa 10-15 minuti |
| ICON-2I, se il riquadro è nella sua griglia | circa 5 minuti |
| ECMWF letto come Shift sulle celle, stabilità, tabella, vento vero, libreria, prova sulla cassaforte, pubblicazione | pochi minuti |

**In tutto circa mezz'ora per area, gratis** nel deposito pubblico.

**Una regola di qualità.** La libreria si pubblica solo se ci sono almeno 70 passaggi del radar e se sulla cassaforte migliora davvero. Altrimenti il flusso lo dice e non pubblica, come oggi per ICON-2I a Venezia e Maiorca.

**Che cosa manca per farlo partire dall'app,** come l'atlante:
- **La lettura del modello dai soli punti di mare deve essere quella di Shift.** Shift decide che cosa è terra con le sue mappe della costa (`europa.rgm`, `world.rgm`, `atlantico.rgm`) e il rilievo con `europa.rgr`. Queste mappe sono nel deposito privato di Shift.
- **Una copia della lettura in Python, con il solo rilievo, non basta.** A Palermo dà gli stessi valori di Shift (scarto massimo 0,003 m/s). In Dalmazia invece, tra le isole, se ne discosta: il 95% delle celle entro 0,46 m/s, qualcuna fino a 3 m/s. È proprio perché Shift usa i contorni della costa e non il rilievo.
- **La strada pulita:** mettere nella release del deposito pubblico `lens` le sole mappe (la costa e il rilievo, fatti da dati pubblici). Poi scrivo in Python, da me, la lettura dei contorni della costa come la fa Shift. Il codice di Shift resta privato.
- **Il resto è lavoro mio e poi di Shift:**
  - il flusso `lens2_area.yml` (riquadro e nome);
  - il controllo che rifaccia uguali le librerie di oggi;
  - in Shift, un tasto "Lens 2 per quest'area" che lo lancia, come per l'atlante.

## 5. File

- **Codice:**
  - `banco/strumenti/vero_tabella.py`: il radar in vento vero sulla tabella;
  - `banco/strumenti/indice_librerie.py`: l'indice dalle librerie;
  - aree nuove in `banco/dati/costa.py`, `stabilita.py`, `tabella_celle.py`, `confronto_metodi.py`.
- **Dati:**
  - release `radar` (le undici aree);
  - `enuove-aperto-*` e `snuove-aperto-*` (ECMWF per Baleari e Venezia);
  - `icon2i-storia`.
