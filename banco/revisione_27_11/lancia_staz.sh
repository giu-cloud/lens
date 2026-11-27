#!/bin/sh
cd /home/claude/lens/banco/strumenti
P="-Dshift.lens.relief=0.25"
for cfg in "_s30:" "_s10:10" "_s5:5"; do
  s=${cfg%%:*}; k=${cfg##*:}
  SUFFISSO=$s PROPS="$P" KM=$k MAX=250 python3 stazioni_c.py palermo corsica-nord cagliari egadi pontine > /home/claude/lavoro/rev/staz$s.a.log 2>&1 &
  SUFFISSO=$s PROPS="$P" KM=$k MAX=250 python3 stazioni_c.py monaco genova malta maddalena-bonifacio > /home/claude/lavoro/rev/staz$s.b.log 2>&1 &
  wait
done
SUFFISSO=_s30 PROPS="$P" MAX=250 python3 stazioni_c.py messina > /home/claude/lavoro/rev/staz_s30.m.log 2>&1 &
SUFFISSO=_v1 PROPS="-Dshift.lens.relief=1" MAX=250 python3 stazioni_c.py palermo corsica-nord monaco messina > /home/claude/lavoro/rev/staz_v1.log 2>&1 &
wait
echo FINITO
