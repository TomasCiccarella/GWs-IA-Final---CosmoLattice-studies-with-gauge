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

## Proposed approach (to confirm before committing compute time)

[CosmoLattice](https://cosmolattice.net/) (Figueroa, Florio, Torrenti &
Valkenburg, JCAP 04 (2021) 035 [arXiv:2006.15122]; code paper
[arXiv:2102.01031]) ships an Abelian-Higgs template out of the box: a
complex scalar charged under a U(1) gauge field, in an expanding FRW
background, with built-in GW spectrum extraction. That makes it the
cleanest lever for this question on a course timescale:

1. Run the **scalar-only control** — the same potential, gauge coupling
   `e → 0`, so the U(1) sector decouples and the run reduces to ordinary
   self-resonance preheating.
2. Run the **Abelian-Higgs case** at a small grid of gauge couplings `e`,
   everything else (potential parameters, lattice size, initial
   conditions, random seed) held fixed.
3. Extract Ω_GW(f) from each run, and compare: peak location and height
   vs. `e`, and whether the gauge-sourced runs leave spectral features
   (a knee, a second bump, a different high-`f` slope) the control does
   not have.

This is a plan, not a result — it still needs: confirming CosmoLattice
builds here (C++17, CMake, MPI, optionally FFTW/HDFI5), picking defensible
potential parameters and a lattice size the machine available can actually
run, and deciding how many gauge-coupling points a course-length budget
affords. Each of those is a choice with a defensible alternative and gets
recorded in `provenance/` when it's made, not narrated here after the fact.

## Repository layout

| | |
|---|---|
| `final-project.html` | the presentable page — what a reader who wasn't in the room sees. Built up as the work progresses, not written at the end. |
| `final-project.pdf` | same content, as a PDF, generated from the HTML once there is something to present |
| `code/` | everything written from scratch: CosmoLattice model/config files, spectrum post-processing, plotting |
| `data/` | run outputs too large or too raw to be a "figure" — pointers/checksums if the actual files don't belong in git |
| `figures/` | every figure that ships in the final page or PDF |
| `papers/` | reference papers (CosmoLattice papers, the preheating-GW literature this builds on) |
| `provenance/` | `claims.yaml`, `numbers.json` — what is asserted and what backs it. Format described in `.claude/provenance/*.md`; created the first time there's a claim, a number, or a figure to record, not before. |
| `.claude/`, `.codex/` | the same provenance gate used in `day5/exercise/` of the course repo: a session-start/stop hook that will not let a turn end with an unrecorded figure or number. |

## Reproducing this

Nothing has been computed yet — this commit is the scaffold, not a result.
Once there's a first run, this section gets replaced with the actual
build/run instructions (compiler, CosmoLattice version/commit, MPI rank
count, how long it took), because a repository that can't be reproduced by
something that has never spoken to me is exactly the failure mode this
course is about.
