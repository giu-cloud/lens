#!/bin/sh
# pezzi separati sopra la base 0 (rilievo 0, ridosso 0.2, costa 0, allungamento 1), e scala del livello con la manopola 0.25
cd /home/claude/lens/banco/cfd
B="-Dshift.lens.relief=0"
run() { # suffisso, props, km
  for a in "messina palermo eolie cagliari" "monaco corsica-nord dalmazia"; do
    SUFFISSO=$1 PROPS="$2" KM=$3 python3 riferimenti.py $a > /home/claude/lavoro/rev/log$1_$(echo $a | cut -d' ' -f1).txt 2>&1 &
  done; wait
}
run _k10 "-Dshift.lens.relief=0.25" 10
run _k15 "-Dshift.lens.relief=0.25" 15
run _pT "$B -Dshift.lens.terrain=1" ""
run _pW "$B -Dshift.lens.wake=0.65" ""
run _pC "$B -Dshift.lens.coast=1" ""
run _pWS "$B -Dshift.lens.wake=0.65 -Dshift.lens.sfacMax=3" ""
run _k5 "-Dshift.lens.relief=0.25" 5
run _k20 "-Dshift.lens.relief=0.25" 20
echo FINITO
