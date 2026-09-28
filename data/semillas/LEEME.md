# Otras semillas del control y de U(1) con vacío transversal (27 de septiembre de 2026)

Para saber cuánto del resultado U(1)/control (bitácora §3.4) es azar. Misma red que las corridas de
la semilla 1234 (N = 64, kIR = 0.5 en unidades del control, VV2); solo cambia `baseSeed`.

| carpeta | modelo | semilla |
|---|---|---|
| `lphi4_N64_kIR0.5_VV2_s2345/`, `..._s3456/` | control | 2345, 3456 |
| `lphi4U1_vacioT_N64_kIR0.5_VV2_s2345/`, `..._s3456/` | U(1) con vacío transversal | 2345, 3456 |

- `correr.sh` genera cada `.in` a partir del de `code/` con otra semilla y las corre en serie
  (`progreso.log` tiene las horas). La semilla 1234 está en `data/convergencia_N/lphi4_N64_kIR0.5_VV2/` y
  `data/lphi4U1_vacioT_N64_kIR0.5_VV2/`.
- Los tiempos de `tiempo.log` de `lphi4U1_..._s2345` (2 h 10 min) no sirven como costo: compartió la
  máquina con compilaciones y pruebas. Sola, una corrida U(1) tarda ~1 h 10 min y una del control ~11 min.
- `resumen.json`: salida de `python3 code/semillas_U1.py` (figura en `figures/semillas/`).
