# Lens 2 contro le stazioni a terra, a tutte le ore (25/11/2026)

Il radar passa solo verso le 05 e le 17 UTC. Per vedere se le correzioni di Lens 2 reggono anche nelle altre ore, e soprattutto le correzioni grandi (Stretto di Messina, Genova, golfo di Palermo), servono misure continue. Le stazioni a terra lo sono, ma molte sono in città, su colline o dietro un capo: leggono il vento del posto, non quello del mare. Qui le stazioni le giudica il radar, settore per settore, e si tengono solo quelle che vanno d'accordo con il mare.

## 1. Le fonti

| Fonte | Che cosa | Stazioni vicine alle aree | Esito |
|---|---|---|---|
| **Météo-France** (data.gouv.fr, licenza aperta) | stazioni orarie, vento medio di 10 minuti a 10 m, controllato | 32 (Corsica, Alpi Marittime) | usata |
| **NOAA GHCNh** (ha preso il posto di ISD nel 2025) | aeroporti e stazioni sinottiche, anche fari e capi dell'Aeronautica | 19 con dati | usata |
| **IEM, Iowa State** | METAR degli aeroporti ogni mezz'ora | 18 | usata |
| **ISPRA, Rete Mareografica Nazionale** | anemometri nei porti (Palermo, Messina, Genova, Napoli, Ponza, Lipari…) | 15 | **inutilizzabile:** i dati aperti finiscono a ottobre 2023, prima del radar |
| **Weather Underground** | stazioni amatoriali | — | **non accessibile:** la chiave gratuita la danno solo a chi ha una stazione propria collegata |
| **Meteonetwork** | stazioni amatoriali | — | **non accessibile per ora:** l'API nuova (v3) chiede una chiave; quelle vecchie non valgono più |

Dopo aver tolto i doppioni (la stessa stazione in due fonti), 57 stazioni candidate.

## 2. Come si scelgono le stazioni buone

1. **Quota:** fuori le stazioni sopra 200 m (colli, passi, fari alti): sentono un altro vento. 22 scartate.
2. **Mare davanti:** servono almeno 3 celle della libreria entro 8 km (5 km nella prova più stretta), lasciando fuori la fascia di 1,5 km dalla costa, dove il radar è spesso vuoto. 14 scartate: aeroporti e paesi nell'interno, come Capodichino, Grazzanise, Decimomannu, Corte.
3. **Il giudizio del radar, per ottante del vento al largo.** Ai passaggi delle settimane di addestramento si confronta la stazione con il radar sulle stesse celle. L'ottante è buono se:
   - ci sono almeno 12 passaggi;
   - la correlazione è almeno 0,6;
   - il rapporto stazione/mare sta tra 0,5 e 1,6, con uno scarto che non supera il 40%.

   Il rapporto di ogni ottante buono è l'**esposizione** della stazione in quel settore. Una stazione di città vale solo per i venti con cui vede il mare. Due stazioni non vanno d'accordo con nessun settore e sono scartate: Cap Sagro e la sinottica di Messina.

**Restano 19 stazioni,** tutte vicine al mare:
- Corsica: Cap Pertusato, Cap Corse, Ile Rousse, Calvi, Bastia;
- Costa Azzurra: Antibes-La Garoupe, Nice;
- aeroporti: Genova, Trapani, Reggio Calabria, Palermo Punta Raisi, Malta;
- capi e isole: Capo Carbonara, Ponza;
- Cicladi: Syros (due stazioni), Mykonos, Paros, Naxos.

Ognuna vale per 1-6 ottanti su 8.

## 3. La prova

- **Quando:** solo le settimane della cassaforte, mai usate per costruire le librerie (librerie di prova senza la cassaforte). Ogni 3 ore, a tutte le ore del giorno.
- **Modelli:**
  - ECMWF dalle corse delle 00 e delle 12, passi da 3 a 12 ore;
  - ICON-2I dalla storia di Open-Meteo, ogni 3 ore, scaricata per questa prova.
- **Il riferimento:** la stazione divisa per la sua esposizione nel settore, cioè il vento sul mare davanti. Il modello si legge sulle stesse celle.
- **Due giudizi:**
  - **assoluto:** errore medio e quadratico, in nodi;
  - **libero:** prima dell'errore quadratico, a ogni stazione e previsione si dà la scala migliore, un solo numero per tutta la serie. Toglie ogni errore costante, anche quello di taratura della stazione, e dice solo se la correzione segue meglio i cambiamenti.
- **Vento al largo:** almeno 4 nodi.
- **Numeri:** 13.765 confronti, di cui 6.262 con ICON-2I alla stessa ora.

## 4. I risultati

### Tutte le stazioni (errore medio / quadratico in nodi, ore con tutti e due i modelli)

| Ore UTC | ECMWF | ECMWF corretto | ICON-2I | ICON-2I corretto | Insieme grezzo | **Insieme corretto** |
|---|---|---|---|---|---|---|
| tutte | −0,79 / 5,28 | −0,48 / 4,91 | −1,53 / 4,96 | −1,24 / 4,78 | −1,09 / 4,82 | **−0,78 / 4,58** |
| notte (00-03) | +0,33 / 4,82 | +0,51 / 4,58 | −0,92 / 4,58 | −0,62 / 4,43 | −0,17 / 4,34 | **+0,06 / 4,23** |
| mattina (06-09) | −0,54 / 5,34 | −0,36 / 4,98 | −1,47 / 4,96 | −1,14 / 4,78 | −0,91 / 4,85 | **−0,67 / 4,63** |
| pomeriggio (12-15) | −2,45 / 5,98 | −2,04 / 5,35 | −2,47 / 5,36 | −2,18 / 5,14 | −2,46 / 5,47 | **−2,10 / 5,02** |
| sera (18-21) | −0,34 / 4,83 | +0,11 / 4,67 | −1,18 / 4,87 | −0,92 / 4,72 | −0,68 / 4,46 | **−0,30 / 4,36** |

Giudizio libero, errore quadratico in nodi:

| Ore | ECMWF → corretto | ICON-2I → corretto | Insieme grezzo → corretto |
|---|---|---|---|
| tutte | 5,03 → 4,85 | 4,66 → 4,64 | 4,58 → 4,47 |
| pomeriggio (12-15) | 5,08 → 4,88 | 4,71 → 4,72 | 4,60 → 4,49 |

Con il raggio più stretto (5 km, 13 stazioni) i numeri vanno nello stesso modo:
- ECMWF da 5,16 a 4,75;
- insieme da 4,61 a 4,36;
- pomeriggio, ECMWF da 5,82 a 5,14.

- **Lens 2 regge a tutte le ore,** anche in quelle che il radar non vede. L'errore quadratico di ECMWF scende del 7%, quello dell'insieme del 5%, e ICON-2I migliora un po'.
- **Nel pomeriggio ECMWF migliora più che nelle altre ore** (−11%). Il pomeriggio è la fascia più lontana dai passaggi del radar.
- **Anche nel giudizio libero si migliora:** non è solo un errore costante tolto.
- **L'insieme corretto resta il migliore a ogni ora.**

### Le correzioni grandi

Errore medio / quadratico in nodi, ECMWF grezzo → corretto:

| Stazione | Fattore medio | Tutte le ore | Pomeriggio (12-15) |
|---|---|---|---|
| **Reggio Calabria (Stretto)** | 1,34 | −4,31 / 6,69 → −0,81 / 5,20 | −7,48 / 8,72 → −3,81 / 6,05 |
| Malta (Gudja) | 1,16 | −1,72 / 4,54 → −0,19 / 4,12 | |
| Mykonos | 1,10 | −3,28 / 5,42 → −1,68 / 4,52 | |
| Trapani | 1,09 | −1,10 / 5,00 → −0,25 / 4,70 | −5,29 / 6,76 → −3,87 / 5,52 |
| Genova | 1,07 | −1,62 / 5,60 → −1,07 / 5,40 | −3,70 / 6,52 → −3,23 / 6,22 |
| Capo Carbonara | 0,88 | +0,70 / 3,67 → −0,43 / 3,52 | peggiora: −0,57 → −1,98 di errore medio |
| Paros | 0,92 | +0,57 / 4,29 → −0,56 / 4,08 | |

Il pomeriggio di Reggio Calabria, Trapani e Genova viene dalla prova a 5 km.

- **Lo Stretto è la conferma più forte.** La correzione toglie quasi tutto l'errore medio a ogni ora, e l'errore quadratico scende di oltre un quinto, nel pomeriggio di quasi un terzo.
- **Dove la correzione alza il vento,** ad aeroporti e isole esposte, la prova le dà ragione dappertutto.
- **Dove lo abbassa,** a Capo Carbonara, Paros e Calvi, il risultato è misto. A Carbonara migliora di notte, la mattina e la sera, ma nel pomeriggio peggiora. A Calvi peggiora (6,63 → 6,71).

### Dove non funziona o non si sa

- **Ponza:** la correzione peggiora un po' (insieme da 3,96 a 4,18). Lo diceva anche la boa di Ponza (sezione 7 di `analoghi_25_11.md`): due prove indipendenti dicono la stessa cosa. Non spengo le Pontine, perché la stazione e la boa dicono qualcosa solo dell'isola, non della costa laziale, ma è un punto da tenere d'occhio.
- **Il golfo di Palermo resta senza prova a tutte le ore.** L'unica stazione utile è Punta Raisi, che vede solo l'angolo ovest della libreria (Isola delle Femmine). Lì la correzione è quasi nulla (fattore 0,99) e neutra. L'anemometro del porto di Palermo (ISPRA) sarebbe quello giusto, ma i suoi dati aperti finiscono nel 2023.

### Il pomeriggio: il limite che resta

Anche corretti, i modelli sottovalutano il vento vicino alla costa nel pomeriggio: −2 nodi in media, fino a −4 nello Stretto e a Trapani, da −2 a −3 a Genova. È la brezza e il riscaldamento del giorno, che il radar delle 05 e delle 17 non vede.

C'è però un'avvertenza onesta. Una stazione a terra di giorno sente anche il rimescolamento sopra la terra calda, quindi una parte di quei 2 nodi può essere della stazione e non del mare. I fari e le isole mostrano un errore del pomeriggio molto più piccolo degli aeroporti, dopo la correzione:
- fari e isole (Cap Pertusato, Cap Corse, Ponza, Capo Carbonara, La Garoupe): da −1,4 a +0,6 nodi;
- aeroporti (Reggio Calabria, Trapani, Calvi, Mykonos, Malta, Genova): da −2,3 a −4,6 nodi.

### La fiducia

La dispersione degli analoghi segue l'errore anche qui, a tutte le ore. Sulle cinque parti uguali, l'errore relativo che resta è:

| Parte | Errore relativo |
|---|---|
| la più sicura | 0,16 |
| seconda | 0,19 |
| terza | 0,24 |
| quarta | 0,28 |
| la meno sicura | 0,34 |

Torna con quanto visto sul radar (da 0,18 a 0,38).

## 5. Che cosa vuol dire

1. **Lens 2 si può passare a Shift così com'è.** Nessuna libreria va spenta.
2. **Le correzioni grandi (Stretto, Malta, Trapani) sono confermate a tutte le ore.** Quelle che abbassano il vento sono meno sicure, e lo sono di più nel pomeriggio.
3. **Il passo naturale dopo è un termine per il pomeriggio,** imparato dalle stazioni dove ci sono, con stagione e differenza di temperatura tra terra e mare. Prima però va separata la parte della stazione da quella del mare, usando i fari e le isole.
4. **Per il golfo di Palermo serve una stazione vera sul mare,** in ordine di bontà:
   - l'anemometro del porto di Palermo, chiedendo a ISPRA i dati dal 2024;
   - le stazioni di Meteonetwork sulla costa (Mondello, Arenella, Sferracavallo), se esposte, con una chiave della v3;
   - una stazione Weather Underground, solo se Giano o un amico ne hanno una collegata.

## 6. File

- **Codice:**
  - `banco/dati/stazioni.py` e il flusso `stazioni.yml`: GHCNh, porti ISPRA, METAR di IEM, Météo-France;
  - `banco/strumenti/stazioni_elenco.py`, `stazioni_prova.py`, `stazioni_sintesi.py`.
- **Dati:** release `stazioni` (serie delle stazioni) e `icon2i-storia` (archivi `*_ogni3ore`).
- **Numeri:** `banco/risultati/stazioni_r8.json` (8 km, 19 stazioni) e `stazioni_r5.json` (5 km, 13 stazioni).
