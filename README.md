# Campos de gauge y el espectro de ondas gravitacionales del recalentamiento

Proyecto final de *Ondas gravitacionales e investigación asistida por IA*
(UBA, 2026). Consigna del curso: [`final-project.html`](https://matiaszaldarriaga.github.io/GW-AI-course/final-project.html).

> **¿No sos de física?** Empezá por [`bitacora.html`](bitacora.html): cuenta
> el proyecto desde cero, qué se hizo, qué se aprendió y por qué importa, con
> un glosario al final.

## En palabras simples

Al terminar la inflación, el campo que la impulsó (el *inflatón*) quedó
oscilando y le pasó su energía a otros campos de forma violenta. Esa etapa se
llama *recalentamiento*, y agitó el espacio lo suficiente como para producir
ondas gravitacionales. Casi todos los estudios suponen que la energía va a
campos "simples" (escalares). Este proyecto simula qué pasa cuando va a
**campos de gauge**, los que transmiten fuerzas como el electromagnetismo, y
mide cuánto cambia la "huella" de ondas gravitacionales que queda.

## La pregunta

El recalentamiento posterior a la inflación genera un fondo estocástico de
ondas gravitacionales, a partir del estrés anisotrópico de configuraciones de
campo clásicas e inhomogéneas en la red. Casi todo lo publicado sobre esto usa
sólo campos escalares. Cuando el contenido de campos incluye un campo de gauge,
porque el inflatón (o un campo espectador) tiene carga bajo algún grupo de
gauge, la energía-momento del propio campo de gauge contribuye a ese estrés. Y
la dinámica de los campos de gauge (tubos de flujo, producción taquiónica de
bosones de gauge, etc.) puede ser muy distinta de la autorresonancia escalar.

**El proyecto pregunta cuánto de la forma del espectro de ondas
gravitacionales (frecuencia y altura del pico, pendientes a cada lado) se
debe al campo de gauge, manteniendo fijo el resto del modelo.**

## Enfoque

[CosmoLattice 2.0](https://cosmolattice.net/) (Baeza-Ballesteros, Figueroa,
Florio, Loayza, Sattler, Torrentí y Urio, arXiv:2607.24978; teoría en
arXiv:2006.15122 y arXiv:2512.15627) trae tres modelos relacionados que
permiten una comparación controlada sin escribir un modelo nuevo:

| modelo | contenido de campos |
|---|---|
| `lphi4` | un escalar real con potencial λφ⁴ y un campo hijo escalar; sin campo de gauge. **El control.** |
| `lphi4U1` | el mismo tipo de potencial, ahora con un escalar complejo cargado bajo un campo de gauge U(1) (electrodinámica escalar / Higgs abeliano). |
| `lphi4SU2U1` | un doblete escalar con carga SU(2)×U(1), parecido al sector electrodébil: dos campos de gauge en lugar de uno. |

Los valores por defecto de los `.in` de CosmoLattice **no** hacen una
comparación justa: los modelos arrancan con energías distintas y con
resonancias de distinta intensidad. El análisis en
[`analisis-parametros.html`](analisis-parametros.html) fija parámetros
comparables, y las decisiones quedaron confirmadas el 22 de septiembre de 2026:

- **λ = 9×10⁻¹⁴** en los tres modelos, con el inflatón arrancando en el estado
  exacto de fin de la inflación.
- **Un único canal resonante con q = 120** en cada modelo: el hijo escalar en
  `lphi4`, el campo U(1) en `lphi4U1`, y los dos campos gauge juntos en
  `lphi4SU2U1` (q_A = q_B = 60, porque según el Art I resuenan con la suma;
  hijos escalares apagados).
- **La misma red física en los tres.** Por un factor √2 en las unidades del
  código, `kIR` y `kCutOff` se multiplican por √2 en los modelos gauge, y `dt`
  y `tMax` se dividen por √2.

Pasos:

1. Correr el **control** (`lphi4`), primero con una corrida piloto de N = 128, y
   validarlo contra el espectro de Dufaux et al. (2007) para q = 120.
2. Correr **`lphi4U1`** y **`lphi4SU2U1`** con los valores emparejados. Chequear
   que la ley de Gauss se cumpla a precisión de máquina desde el inicio y que
   el inflatón oscile igual que en el control.
3. Extraer Ω_GW(f) de cada corrida y comparar: posición y altura del pico, y si
   los modelos con gauge dejan rasgos que el control no tiene (un quiebre, un
   segundo pico, otra pendiente a alta frecuencia).

Hay un intento previo, sin terminar, con este mismo modelo en mi trabajo de
tesis. No se reutiliza acá directamente: un repositorio que otra persona clona
desde GitHub tiene que ser reproducible por sí solo, sin depender de archivos
externos. Pero conviene comparar contra él cuando haya un resultado.

## Organización del repositorio

| | |
|---|---|
| `bitacora.html` | **para empezar si no sos de física**: diario de trabajo en lenguaje llano, con lo hecho, lo aprendido y un glosario. Se actualiza a medida que avanza el trabajo. |
| `final-project.html` | la página de presentación, para un lector que no estuvo en clase. Se arma a medida que avanza el trabajo, no al final. |
| `what-is-cosmolattice.html` | introducción a CosmoLattice para quien nunca lo usó: qué es, cómo pone campos escalares y de gauge U(1)/SU(2) en una red, y las ecuaciones que generan las ondas gravitacionales, cada una citada con su número de ecuación. |
| `bases-teoricas-modelos-gauge.html` | un nivel más abajo: variables de programa, condiciones iniciales y ley de Gauss discreta en los tres modelos, cada paso verificado por `code/verificar_variables_de_programa.py`. |
| `analisis-parametros.html` | qué parámetros hacen comparables a los tres modelos: los valores que usan los autores y por qué, λ según las observaciones, un análisis de Floquet de q y la tabla recomendada. Los números salen de `code/analisis_parametros.py`. |
| `final-project.pdf` | el mismo contenido en PDF; se genera desde el HTML cuando haya algo para presentar. |
| `code/` | todo lo escrito desde cero: verificaciones, análisis, archivos de configuración de las corridas, post-procesamiento y gráficos. |
| `data/` | salidas de las corridas demasiado grandes o crudas para ser una figura (punteros y checksums si los archivos no van en git). |
| `figures/` | todas las figuras que aparecen en las páginas o en el PDF. |
| `bibliografía/` | los papers de referencia: los de CosmoLattice (código, teoría, GWs) y la literatura sobre campos de gauge en el recalentamiento. `bibliografía/BIBLIOGRAPHY.md` explica qué es cada uno y para qué está. |
| `CosmoLattice/` | el código de CosmoLattice, clonado de upstream. No se commitea (ver "Cómo reproducirlo"): se compila desde la fuente cada vez. |
| `provenance/` | `claims.yaml` y `numbers.json`: qué se afirma y qué lo respalda. El formato está en `.claude/provenance/*.md`. |
| `.claude/`, `.codex/` | el mismo control de procedencia que en `day5/exercise/` del repo del curso: un hook de inicio y fin de sesión que no deja terminar un turno con una figura o un número sin registrar. |

## Cómo reproducirlo

Todavía no se corrió ninguna simulación física; es el próximo paso. Lo que
está hecho es el entorno: CosmoLattice compila sin problemas acá, en tres de
sus modelos.

```bash
git clone https://github.com/cosmolattice/cosmolattice.git CosmoLattice
cd CosmoLattice
mkdir build_lphi4 && cd build_lphi4
cmake -DMODEL=lphi4 -DOPENMP=ON ..
make cosmolattice -j"$(nproc)"
# repetir con -DMODEL=lphi4U1 y -DMODEL=lphi4SU2U1, cada uno en su propio directorio de build
```

- Commit de upstream: `acc8278d8832890754a1df16aec9eab5e1867c5c` (2026-08-04),
  `https://github.com/cosmolattice/cosmolattice`. Ya incluye el fix
  `f6b9c267` del medidor de la ley de Gauss en SU(2). El `FetchContent` de
  CMake descarga automáticamente al configurar el backend TempLat (fijado en
  `v1.0.2` por el propio `CMakeLists.txt` de CosmoLattice) y Kokkos. Hace falta
  conexión a internet una vez; no se incluyen en el repo.
- Herramientas usadas: `g++` 13.3.0 (Ubuntu 24.04), CMake 4.0.3, GNU Make 4.3,
  compilado con `-DOPENMP=ON` (Kokkos detectó OpenMP como backend de CPU; sin
  GPU ni MPI).
- Los tres builds (`lphi4`, `lphi4U1`, `lphi4SU2U1`) compilaron a un binario
  funcional con una misma advertencia inofensiva (`-Wshadow` sobre `FloatType`
  en `abstractmodel.h`, presente en los tres): un choque de nombres dentro de
  la jerarquía de templates de CosmoLattice, no algo introducido acá.
- Prueba de humo de `lphi4`: `./lphi4 input=../../models/parameter-files/lphi4.in N=16 tMax=0.5`
  terminó bien (código de salida 0) y escribió los archivos esperados
  (`average_*.txt`, `spectra_*.txt`). Es un chequeo del build, no un
  resultado físico: `N=16` y `tMax=0.5` son demasiado chicos para significar
  algo, y no se guardó nada de esa corrida.
- Las verificaciones y el análisis se reproducen con
  `python3 code/verificar_variables_de_programa.py` (necesita `sympy`) y
  `python3 code/analisis_parametros.py [--tabla | --figura]` (necesita
  `numpy`, `scipy` y `matplotlib`).
- El código de CosmoLattice y los directorios de build se clonan de cero y no
  se commitean (`.gitignore` excluye `CosmoLattice/`): es una dependencia
  externa, fijada por el hash de commit de arriba, no algo que escribió este
  proyecto.

Cuando haya una primera corrida física, esta sección va a incluir los
parámetros usados, el tiempo de cómputo y cómo se verificaron las salidas.
