# SU(2)×U(1) con vacío transversal (27-28 de septiembre de 2026)

La corrida larga de `lphi4SU2U1`: N = 64, kIR = 0.5 en unidades del control (0.707107 en unidades
gauge), VV2, semilla 1234, q_A = q_B = 60, con `ICtype_U1` e `ICtype_SU2 = RandomWithMatterTransverseVacuum`
(parches `code/parches/u1_` y `su2_vacio_transversal.patch`, build `CosmoLattice/build_lphi4SU2U1_tv`).

- `correr.sh` esperó a que terminaran las semillas (`data/semillas/`) y la largó: 23:18 del 27/9 a 08:40
  del 28/9, 9 h 22 min con 8 núcleos (1,12 s/paso), 249 MB, exit 0 (`progreso.log`, `tiempo.log`).
- `salida.log` muestra la matriz de masas y los autoestados que usó el parche (fotón, Z, W, W).
- `analisis.json`: salida de `python3 code/analisis_SU2U1.py` (figuras en
  `figures/lphi4SU2U1_vacioT_N64_kIR0.5_VV2/`).
- `crecimiento_por_modo.json`: salida de `python3 code/crecimiento_por_modo_U1.py lphi4SU2U1_vacioT_N64_kIR0.5_VV2`
  (μ por modo del espectro magnético del U(1), que en este modelo es la mitad "Z + fotón").
- Ojo: el espectro `spectra_norm_SU2_0.txt` es el de |B| y |E| (escalares), no el de cada modo de B^a.
  La energía magnética de SU(2) tiene un piso de redondeo (~10⁻¹² de la energía total) mientras los
  campos son de vacío. Ver bitácora §3.5 y §3.7.
