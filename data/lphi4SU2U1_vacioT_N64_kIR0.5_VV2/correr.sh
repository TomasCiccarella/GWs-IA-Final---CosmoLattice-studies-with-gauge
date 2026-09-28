#!/bin/bash
# Corrida larga de SU(2)×U(1) con vacío transversal (parches u1_ y su2_vacio_transversal, build_lphi4SU2U1_tv).
# Espera a que terminen las semillas (data/semillas/correr.sh, PID 19584) para tener los 8 núcleos libres.
cd "$(dirname "$0")"
while kill -0 19584 2>/dev/null; do sleep 60; done
echo "arranque: $(date)" > progreso.log
/usr/bin/time -v -o tiempo.log /home/tomy/Desktop/UBA/Curso_GWsIA/GW-AI-course/Final/CosmoLattice/build_lphi4SU2U1_tv/lphi4SU2U1 input=lphi4SU2U1_vacioT_N64_kIR0.5_VV2.in > salida.log 2>&1
echo "terminada (exit $?): $(date)" >> progreso.log
