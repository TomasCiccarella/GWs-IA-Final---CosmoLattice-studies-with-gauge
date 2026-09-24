# Parches a CosmoLattice

`CosmoLattice/` no está versionado en este repo (es un clon de upstream, commit `acc8278d`).
Los cambios que le hacemos quedan acá como parches.

## `u1_vacio_transversal.patch`

Agrega la condición inicial `ICtype_U1 = RandomWithMatterTransverseVacuum` (o `7`) para modelos
con un campo U(1) acoplado a escalares complejos (`lphi4U1`).

- **Por qué:** la condición por defecto (`RandomWithMatter`, Art I §7.2 y paper del código
  2102.01031 ecs. 97–102) arranca con A = 0 y solo el campo eléctrico *longitudinal* que fija la
  ley de Gauss. Los modos *transversales*, que son los que resuenan y emiten GWs, arrancan en cero
  y solo se siembran a segundo orden. El campo hijo χ del control, en cambio, arranca con
  fluctuaciones de vacío. Ver `bitacora.html#condicion-inicial-gauge`.
- **Qué hace:** lo mismo que `RandomWithMatter` (escalares + Gauss, con las mismas semillas) y
  además suma fluctuaciones de vacío a A y E en las dos polarizaciones transversales, con la
  función `planeWaves` que ya trae CosmoLattice. Son transversales respecto del momento de red
  hacia atrás, así que no cambian la divergencia de red y la ley de Gauss sigue igual.
  - Masa efectiva de cada polarización: m² = 2 (gQ)² |Φ₀|² / ω★² (unidades de programa, a = 1),
    igual que χ se inicializa con su masa efectiva.
  - `aDot = 0` en el momento: el campo gauge es conforme, no lleva el término −a′/a·A de los escalares.
  - `planeWaves` ganó un parámetro opcional `mass2` (por defecto 0, como antes).
  - `extrafields.h`: la memoria auxiliar también se reserva cuando la condición la pide el `.in`.
- **Aplicar:** `cd CosmoLattice && git apply ../code/parches/u1_vacio_transversal.patch`,
  y compilar en `build_lphi4U1_tv` (`cmake .. -DMODEL=lphi4U1 ...`).
