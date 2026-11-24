#!/bin/sh
# Il metro di misura del laboratorio: ECMWF 0,25 letto come Shift (normale e dai soli punti di mare) contro il radar,
# tutte le aree. Rifà la tabella di routing/atlante_costa_risultati.md (24/11).
#   DATI=<cartella dei costa_*.npz> ECMWF=<cartella degli ecmwf_*.grb2> USCITA=<cartella> CORE=<core di Shift .jar> sh base_ecmwf.sh
# Servono: kotlinc 2.0.21 (in /tmp/kotlinc), il nucleo di Shift compilato (CORE) e le mappe di Shift in
# /home/claude/Routing/app/src/main/assets/maps (il percorso è scritto in ShiftSar.kt).
set -e
cd "$(dirname "$0")/../costa"
[ -f shiftsar.jar ] || /tmp/kotlinc/bin/kotlinc ShiftSar.kt -cp "$CORE" -d shiftsar.jar
UT="$CORE" NOLENS=1 python3 costa_ecmwf.py "$DATI" - "$ECMWF" "$USCITA" egadi eolie palermo messina malta napoli cagliari maddalena-bonifacio toscano genova
python3 costa_pool.py "$USCITA" cmems shift mare
