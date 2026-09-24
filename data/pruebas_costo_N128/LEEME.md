# Pruebas de costo a N = 128 (23 de septiembre de 2026)

Corridas cortas (100–200 pasos) para medir cuánto tardan y cuánta memoria usan los tres
modelos con N = 128, GWs activadas y el integrador VV2. **No son resultados físicos**: sólo
miden costo. Los parámetros son los de la tabla de `analisis-parametros.html#tabla`
(traducidos por √2 en los modelos gauge), con `tMax` corto.

| modelo | segundos por paso | corrida completa (~30 000 pasos) | memoria máxima |
|---|---|---|---|
| `lphi4` (control, VV2) | 0,16 | ~1 h 20 min | 318 MB |
| `lphi4U1` | 1,04 | ~8,5 h | 543 MB |
| `lphi4SU2U1` | 10,4 | ~86 h (~3,5 días) | 1,56 GB |

- Segundos por paso: diferencia de las marcas de tiempo de `salida.log` entre el paso 0 y el
  paso 100. Memoria y uso de CPU (~790 %, los 8 núcleos con OpenMP): `tiempo.log` (`/usr/bin/time -v`).
- **Por qué VV2 y no LF:** CosmoLattice no permite LF con campos U(1) y GWs a la vez
  (`include/CosmoInterface/evolvers/evolver.h`, línea 79): el kernel de las GWs necesita los
  campos eléctrico y magnético sincronizados. En el control VV2 cuesta lo mismo que LF.
- **Ley de Gauss en t = 0:** violación relativa ~2×10⁻¹¹ en U(1) y ~4×10⁻¹⁰ en SU(2)
  (`average_gauss_*.txt`).
- La primera prueba con LF abortó con `InvalidEvolverTypeGW`; de ahí salió lo anterior.
