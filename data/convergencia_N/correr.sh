#!/bin/bash
# Corre en serie las corridas de convergencia (cada una usa los 8 núcleos con OpenMP).
BIN=/home/tomy/Desktop/UBA/Curso_GWsIA/GW-AI-course/Final/CosmoLattice/build_lphi4/lphi4
cd "$(dirname "$0")"
for n in lphi4_N32_VV2 lphi4_N48_VV2 lphi4_N64_VV2 lphi4_N64_kIR0.5_VV2 lphi4_N96_VV2; do
  (cd $n && /usr/bin/time -v -o tiempo.log $BIN input=$n.in > salida.log 2>&1)
  echo "$n terminada: $(date)"
done
