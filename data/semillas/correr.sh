#!/bin/bash
# Otras semillas para U(1) con vacío transversal y su control (N = 64, kIR = 0.5, VV2).
# Cada .in es el de code/ con otra baseSeed; la semilla 1234 ya está en
# data/lphi4U1_vacioT_N64_kIR0.5_VV2/ y data/convergencia_N/lphi4_N64_kIR0.5_VV2/.
# Corren en serie: cada una usa los 8 núcleos con OpenMP.
RAIZ=/home/tomy/Desktop/UBA/Curso_GWsIA/GW-AI-course/Final
cd "$(dirname "$0")"
for s in 2345 3456; do
  for m in lphi4:lphi4_N64_kIR0.5_VV2:build_lphi4/lphi4 lphi4U1:lphi4U1_vacioT_N64_kIR0.5_VV2:build_lphi4U1_tv/lphi4U1; do
    IFS=: read -r modelo base bin <<< "$m"
    n=${base}_s$s
    mkdir -p $n
    sed "s/^baseSeed = .*/baseSeed = $s/" $RAIZ/code/$base.in > $n/$n.in
    (cd $n && /usr/bin/time -v -o tiempo.log $RAIZ/CosmoLattice/$bin input=$n.in > salida.log 2>&1)
    echo "$n terminada: $(date)"
  done
done
