#!/bin/bash
# Varianti della lettura "solo mare" di ECMWF, senza Lens, su tutte le aree.
cd "$(dirname "$0")"
AREE="palermo egadi eolie napoli toscano maddalena-bonifacio messina malta cagliari genova"
run() { mkdir -p outM/$1; NOLENS=1 UT=$PWD/../ut2 SEA_PROPS="$2" python3 costa_ecmwf.py dati atlanti ecmwf outM/$1 $AREE > outM/$1/log.txt 2>&1; echo "$1 fatto"; }
run point0 '{"sea.mask":"point","sea.landw":0}'
run frac50 '{"sea.mask":"frac","sea.frac":0.5,"sea.landw":0}'
run point25 '{"sea.mask":"point","sea.landw":0.25}'
run frac30 '{"sea.mask":"frac","sea.frac":0.3,"sea.landw":0}'
run point50 '{"sea.mask":"point","sea.landw":0.5}'
run frac70 '{"sea.mask":"frac","sea.frac":0.7,"sea.landw":0}'
