#!/bin/bash
# Otra semilla para SU(2)×U(1) con vacío transversal (N = 64, kIR = 0.5, VV2, build_lphi4SU2U1_tv).
# La semilla 1234 ya está en data/lphi4SU2U1_vacioT_N64_kIR0.5_VV2/. Usa los 8 núcleos (~9 h 20 min).
RAIZ=/home/tomy/Desktop/UBA/Curso_GWsIA/GW-AI-course/Final
cd "$(dirname "$0")"
base=lphi4SU2U1_vacioT_N64_kIR0.5_VV2
for s in 2345; do
  n=${base}_s$s
  mkdir -p $n
  sed "s/^baseSeed = .*/baseSeed = $s/" $RAIZ/code/$base.in > $n/$n.in
  echo "$n arranque: $(date)" >> progreso.log
  (cd $n && /usr/bin/time -v -o tiempo.log $RAIZ/CosmoLattice/build_lphi4SU2U1_tv/lphi4SU2U1 input=$n.in > salida.log 2>&1)
  echo "$n terminada (exit $?): $(date)" >> progreso.log
done
