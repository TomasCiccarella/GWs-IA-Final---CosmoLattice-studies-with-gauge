# Pruebas del parche de vacío transversal para SU(2) (27 de septiembre de 2026)

Corridas cortas de `lphi4SU2U1` (N = 64, kIR = 0.5 en unidades del control, misma semilla 1234),
hechas para validar `code/parches/su2_vacio_transversal.patch` antes de la corrida larga.
**No son resultados físicos.** Explicación en `bitacora.html#vacio-su2`.

| carpeta | qué es | hasta t (programa) |
|---|---|---|
| `defecto/` | condición inicial de CosmoLattice (links = 1, solo E longitudinal) | 1,06 |
| `vacioT/` | vacío transversal en la base de autoestados de masa (versión final del parche) | 0,71 |
| `vacioT_amp1000/` | lo mismo con las fluctuaciones transversales ×1000 (`COSMOLATTICE_VACIO_T_AMP=1000`) | 0,71 |
| `v1_sin_mezcla/` | primera versión del parche: A y cada color con su masa diagonal (120), sin la mezcla | 1,06 |
| `U1ref_defecto/`, `U1ref_vacioT/` | referencia: `lphi4U1` (parche U(1) ya validado), misma red y salida densa | 0,71 |

- `validacion.json`: salida de `python3 code/validar_vacio_su2.py` (compara `defecto/` con `vacioT/`).
- Las pruebas corrieron con 4 hilos mientras corrían las semillas (`data/semillas/`), así que los
  tiempos de `tiempo.log` no sirven como medida de costo.
- `defecto/` y `v1_sin_mezcla/vacioT/` se cortaron a mano en t ≈ 1,06 (por eso no tienen `tiempo.log`
  completo); los datos hasta ahí están bien.

Resultado (ver `provenance/numbers.json`, claves `SU2U1_*` y `SU2_vT_*`):

- La matriz de masas mezcla A con el color 1: fotón (A + B¹)/√2 sin masa, Z (A − B¹)/√2 con m² = 240,
  W (B², B³) con m² = 120.
- Energía agregada en t = 0: magnética 1,03 de la analítica, eléctrica 0,90 (igual en U(1) y SU(2)).
- Con la mezcla, E·a⁴ y B·a⁴ quedan casi constantes, como en la referencia `lphi4U1`. Sin la mezcla
  (`v1_sin_mezcla/`), B·a⁴ del U(1) subía ×3,6: el vacío no correspondía a las masas reales.
- La energía magnética de SU(2) se calcula con 2 − tr(plaqueta), que con campos de vacío pierde
  precisión (redondeo). Con amplitud física B_SU2·a⁴ parece subir a 1,35; con ×1000 queda en 0,86–1,0.
- Ley de Gauss SU(2): 3,5×10⁻⁷ con vacío transversal (1,5×10⁻⁹ sin él); escala como amplitud².
