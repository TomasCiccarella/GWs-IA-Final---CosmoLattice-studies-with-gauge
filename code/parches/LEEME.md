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
  fluctuaciones de vacío. Ver `paginas/bitacora.html#condicion-inicial-gauge`.
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

## `su2_vacio_transversal.patch`

Va **después** de `u1_vacio_transversal.patch` (lo extiende). Agrega la opción
`ICtype_SU2 = RandomWithMatterTransverseVacuum` (o `2`; por defecto `RandomWithMatter`, igual que
antes) para modelos con un campo SU(2) acoplado a un doblete (`lphi4SU2U1`).

- **Por qué:** SU(2) arranca con el mismo problema que U(1). En `su2initializer.h`, CosmoLattice pone
  los links en la identidad (B^a = 0 exacto) y solo el campo eléctrico longitudinal que fija la ley de
  Gauss. Los modos transversales, los que resuenan y emiten GWs, arrancan en cero. Ver
  `paginas/bitacora.html#vacio-su2`.
- **La mezcla (tipo Z/fotón):** con el doblete en el fondo, la matriz de masas de los campos gauge
  mezcla el U(1) A con el color de SU(2) que apunta en n^a ∝ Φ†σ^aΦ. Si también se pide
  `ICtype_U1 = RandomWithMatterTransverseVacuum`, el parche genera el vacío en la base de
  **autoestados de masa**. Con q_A = q_B = 60 son:
  - un "fotón" (A + B¹)/√2 **sin masa**,
  - una "Z" (A − B¹)/√2 con m² = 240,
  - dos "W" (B², B³) con m² = 120 (unidades de programa).
  Después `U1Initializer` conserva esa parte transversal de A (no genera otra).
  Sin `ICtype_U1` transversal se incluyen solo los tres colores de SU(2).
- **La matriz de masas no está escrita a mano:** se calcula al arrancar a partir del doblete
  homogéneo, con las mismas primitivas que usan `MatterCurrents::U1Current` y `SU2Current`
  (`scalar_prod`, `i_sigma`, `complexPhase`) linealizadas en los campos gauge. Así sigue las
  convenciones del código para el doblete, las matrices de Pauli y los signos de los links. El
  programa verifica que salga simétrica (si no, aborta), la diagonaliza (Jacobi) e imprime en
  `salida.log` la matriz, las masas y las componentes de cada autoestado.
- **Convenciones de SU(2)** (de `FieldFunctionals::B2SU2`, `pi2SU2`, `SU2Kernels` y el mapa exponencial
  del integrador): las componentes del álgebra del link valen g·dx·B^b/2 y `piSU2` guarda E^b/2.
  `unitarize()` recalcula c₀ = √(1 − |c|²).
- **Semillas:** el ruido de CosmoLattice se siembra con `baseSeed` + el nombre del campo auxiliar. Cada
  autoestado usa auxiliares con nombre propio (`SU2ICvacModo0…3`), así el ruido es independiente entre
  ellos y del de los escalares.
- **Ley de Gauss:** la parte transversal tiene divergencia de red nula, pero el término no abeliano
  g[B, E] hace que la ley de Gauss de SU(2) cambie a segundo orden en estas fluctuaciones: la violación
  relativa pasa de 1,5×10⁻⁹ a 3,5×10⁻⁷ (N = 64). La de U(1) no cambia (10⁻¹⁰).
- **Qué cambia en U(1):** `ICtype_U1` transversal en un modelo con doblete deja la parte transversal a
  cargo de `SU2Initializer` cuando también se pide `ICtype_SU2` transversal. Si se pide solo la de U(1),
  la masa suma los dobletes (m_A² = 2 (g_A Q_A)² |Φ₀|² / ω★²) pero ignora la mezcla.
- **Gancho solo para pruebas:** la variable de entorno `COSMOLATTICE_VACIO_T_AMP` multiplica las
  fluctuaciones transversales (U(1) y SU(2)); sin definirla vale 1 y no cambia nada. Sirvió para
  comprobar la dinámica por encima del piso de redondeo del diagnóstico de energía magnética de SU(2),
  que usa 2 − tr(plaqueta) y con campos de vacío resta números de orden 1 que difieren en ~10⁻¹⁴.
- **Validación:** `code/validar_vacio_su2.py` sobre `data/prueba_su2_vacioT/` (ver su `LEEME.md`).
- **Aplicar:** `cd CosmoLattice && git apply ../code/parches/u1_vacio_transversal.patch && git apply
  ../code/parches/su2_vacio_transversal.patch`, y compilar en `build_lphi4SU2U1_tv`
  (`cmake .. -DMODEL=lphi4SU2U1 -DOPENMP=ON && make lphi4SU2U1`).
