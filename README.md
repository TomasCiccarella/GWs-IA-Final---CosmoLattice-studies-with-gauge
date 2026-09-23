# Gauge fields and the gravitational-wave spectrum from preheating

Final project for *Gravitational Waves and AI-Assisted Research* (UBA,
2026). Course brief: [`final-project.html`](https://matiaszaldarriaga.github.io/GW-AI-course/final-project.html).

## The question

Post-inflationary preheating sources a stochastic gravitational-wave
background from the anisotropic stress of classical, inhomogeneous field
configurations on the lattice. Almost everything published on this uses
scalar fields alone. When the field content includes a gauge field —
because the inflaton (or a spectator) is charged under some gauge group —
the gauge field's own energy-momentum contributes to that stress, and
gauge-field dynamics (flux tubes, (p)reheating via tachyonic gauge-boson
production, etc.) can differ sharply from scalar self-resonance.

**This project asks how much of the GW spectrum's shape — peak frequency,
peak amplitude, spectral slopes either side of the peak — is attributable
to the gauge field, holding the rest of the model fixed.**

## Proposed approach (physics choices below still to confirm)

[CosmoLattice 2.0](https://cosmolattice.net/) (Baeza-Ballesteros, Figueroa,
Florio, Loayza, Sattler, Torrentí & Urio, arXiv:2607.24978; theory in
arXiv:2006.15122 and arXiv:2512.15627) ships three related built-in models
that make a controlled comparison possible without writing a new model
from scratch:

| model | field content |
|---|---|
| `lphi4` | a real scalar with a λφ⁴ potential — no gauge field. The control. |
| `lphi4U1` | the same kind of potential, now a complex scalar charged under a U(1) gauge field (scalar electrodynamics / Abelian Higgs). |
| `lphi4SU2U1` | an SU(2)×U(1)-charged scalar doublet (electroweak-like) — two gauge fields instead of one. |

1. Run the **scalar-only control** (`lphi4`).
2. Run **`lphi4U1`** at a small grid of U(1) gauge couplings `gU1s`,
   everything else (potential parameters, lattice size, initial
   conditions, random seed) held fixed.
3. If time allows, run **`lphi4SU2U1`** as a second gauge case — a
   qualitatively different gauge sector (non-Abelian) rather than just a
   different coupling strength. An earlier, unfinished attempt at this
   same model exists in my thesis work; it isn't reused here directly
   (a repository someone else clones from GitHub has to be reproducible
   on its own, not depend on files that live outside it), but it's worth
   checking any result here against once there is one.
4. Extract Ω_GW(f) from each run, and compare: peak location and height,
   and whether the gauge-sourced runs leave spectral features (a knee, a
   second bump, a different high-`f` slope) the control does not have.

Still open, and each is a choice that gets recorded in `provenance/` when
it's actually made, not narrated here in advance: the potential parameters
(`lambda`, `q`, initial amplitudes/momenta) and lattice size (`N`, `kIR`,
`dt`, `tMax`) — the defaults in each model's `.in` file are tuned for
*some* physics, not necessarily one that shows a clean gauge-field effect
— and how many gauge-coupling points a course-length compute budget
affords.

## Repository layout

| | |
|---|---|
| `bitacora.html` | **start here if you're not a physicist** — a plain-language work log in Spanish: what the project asks, what was done, what was learned and why it matters, with a glossary. Updated as the work progresses. |
| `final-project.html` | the presentable page — what a reader who wasn't in the room sees. Built up as the work progresses, not written at the end. |
| `what-is-cosmolattice.html` | a primer for a cosmologist who's never used CosmoLattice: what it is, how it puts scalar and U(1)/SU(2) gauge fields on a lattice, and the GW-sourcing equations, every one cited to an equation number in `bibliografía/`. |
| `analisis-parametros.html` | (Spanish) which parameters make `lphi4`, `lphi4U1` and `lphi4SU2U1` comparable: every value the authors use and why, observational λ, a Floquet analysis of q, and the recommended table. Numbers from `code/analisis_parametros.py`. |
| `bases-teoricas-modelos-gauge.html` | the next level down, in Spanish: program variables, initial conditions and the discrete Gauss law in `lphi4`/`lphi4U1`/`lphi4SU2U1`, each step checked symbolically by `code/verificar_variables_de_programa.py`. |
| `final-project.pdf` | same content, as a PDF, generated from the HTML once there is something to present |
| `code/` | everything written from scratch: CosmoLattice model/config files, spectrum post-processing, plotting |
| `data/` | run outputs too large or too raw to be a "figure" — pointers/checksums if the actual files don't belong in git |
| `figures/` | every figure that ships in the final page or PDF |
| `bibliografía/` | reference papers — the CosmoLattice code/theory/GW papers and the gauge-field-at-preheating literature this builds on. See `bibliografía/BIBLIOGRAPHY.md` for what each one is and why it's here. |
| `CosmoLattice/` | the code itself, `git clone`d from upstream — not committed (see "Reproducing this"), rebuilt from source every time. |
| `provenance/` | `claims.yaml`, `numbers.json` — what is asserted and what backs it. Format described in `.claude/provenance/*.md`; created the first time there's a claim, a number, or a figure to record, not before. |
| `.claude/`, `.codex/` | the same provenance gate used in `day5/exercise/` of the course repo: a session-start/stop hook that will not let a turn end with an unrecorded figure or number. |

## Reproducing this

No simulation has run yet — that's the next step. What's done so far is
the environment: CosmoLattice builds cleanly here, in three of its
built-in models.

```bash
git clone https://github.com/cosmolattice/cosmolattice.git CosmoLattice
cd CosmoLattice
mkdir build_lphi4 && cd build_lphi4
cmake -DMODEL=lphi4 -DOPENMP=ON ..
make cosmolattice -j"$(nproc)"
# repeat with -DMODEL=lphi4U1 and -DMODEL=lphi4SU2U1 in their own build dirs
```

- Upstream commit: `acc8278d8832890754a1df16aec9eab5e1867c5c` (2026-08-04),
  `https://github.com/cosmolattice/cosmolattice`. CMake's `FetchContent`
  pulls the TempLat backend (pinned by CosmoLattice's own `CMakeLists.txt`
  to `v1.0.2`) and Kokkos automatically at configure time — both need
  network access once, not vendored here.
- Toolchain used: `g++` 13.3.0 (Ubuntu 24.04), CMake 4.0.3, GNU Make 4.3,
  built with `-DOPENMP=ON` (Kokkos auto-detected OpenMP as the CPU
  backend; no GPU, no MPI).
- All three builds (`lphi4`, `lphi4U1`, `lphi4SU2U1`) compiled to a working
  binary with one identical, harmless warning (`-Wshadow` on `FloatType`
  in `abstractmodel.h`, present in all three — a naming collision inside
  CosmoLattice's own template hierarchy, not something introduced here).
- `lphi4` was smoke-tested: `./lphi4 input=../../models/parameter-files/lphi4.in N=16 tMax=0.5`
  ran to completion (exit 0) and wrote the expected `average_*.txt` /
  `spectra_*.txt` output files. This is a build check, not a physics
  result — `N=16` and `tMax=0.5` are far too small/short to mean anything,
  and nothing from this run is kept.
- CosmoLattice's own source and every build directory are `git clone`d
  fresh, not committed (`.gitignore` excludes `CosmoLattice/`): it's an
  external dependency pinned by commit hash above, not part of what this
  project wrote.

Once there's a first physics run, this section grows to cover the actual
run parameters, wall-clock time, and how output was checked.
