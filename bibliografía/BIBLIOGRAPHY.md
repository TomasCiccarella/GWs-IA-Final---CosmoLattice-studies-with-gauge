# Bibliography

Every PDF here was pulled from arXiv on 2026-09-22, sourced from
[cosmolattice.net](https://cosmolattice.net) (its documentation, GW manual
page, gauge-field manual page, and citation page) plus a literature search
for prior work specifically on gauge fields and gravitational waves at
preheating. Grouped by what each paper is for in this project, not by
publication date.

## CosmoLattice itself — code, backend, and how to cite it

| file | reference | why it's here |
|---|---|---|
| `figueroa2021_cosmolattice_code_2102.01031.pdf` | Figueroa, Florio, Torrentí, Valkenburg, *CosmoLattice — User Manual*, arXiv:2102.01031 | The v1.0 code paper / user manual: installation, input files, how to define a model. The base reference for everything the build in `../CosmoLattice/` does. |
| `baezaballesteros2026_cosmolattice2.0_2607.24978.pdf` | Baeza-Ballesteros, Figueroa, Florio, Loayza, Sattler, Torrentí, Urio, *CosmoLattice 2.0*, arXiv:2607.24978 | What's actually cloned and built here — v2.0 restructured the repo (`models/` at top level, C++20, TempLat backend) and this is its own code paper. |
| `florio2026_templat_2607.24908.pdf` | Florio & Sattler, *TempLat: a versatile C++ engine for lattice field theories*, arXiv:2607.24908 | TempLat is the backend CosmoLattice 2.0 links against (pulled automatically by CMake's `FetchContent`). Needed to actually read the build, not just run it. |
| `figueroa2021_art1_theory_lattice_2006.15122.pdf` | Figueroa, Florio, Torrentí, Valkenburg, *The Art of Simulating the Early Universe — Part I: Integration Techniques and Canonical Cases*, arXiv:2006.15122 | The theory/lattice-discretization dissertation behind the code: how scalar **and gauge** fields are put on a lattice, the integrators, the canonical example models (this is where `lphi4U1`/`lphi4SU2U1`'s lattice equations come from). |
| `baezaballesteros2026_art2_noncanonical_gw_2512.15627.pdf` | Baeza-Ballesteros, Figueroa, Florio, Lizarraga, Loayza, Marschall, Opferkuch, Stefanek, Torrentí, Urio, *The Art of Simulating the Early Universe — Part II: Non-canonical Cases & Gravitational Waves*, JCAP 06:087 (2026), arXiv:2512.15627 | **The GW methodology paper.** Documents the current (5-dof) algorithm CosmoLattice 2.0 uses to evolve the transverse-traceless GW tensor and extract Ω_GW(f), including the anisotropic-stress contributions from gauge fields. This is what any GW number out of this project has to cite as `produced_by`/method. |

## Gauge fields and GW at preheating — the physics literature this project sits in

| file | reference | why it's here |
|---|---|---|
| `garciabellido2007_stochastic_gw_hybrid_preheating_astroph0701014.pdf` | García-Bellido & Figueroa, *A stochastic background of gravitational waves from hybrid preheating*, arXiv:astro-ph/0701014 | Earliest of this group's papers proposing hybrid-inflation (gauge-charged waterfall field) preheating as a GW source — the scenario this project's control/gauge comparison descends from. |
| `garciabellido2008_gw_background_hybrid_inflation_0707.0839.pdf` | García-Bellido, Figueroa & Sastre, *A Gravitational Wave Background from Reheating after Hybrid Inflation*, Phys. Rev. D 77, 043517 (2008), arXiv:0707.0839 | Follow-up with the full three-stage picture (tachyonic preheating → bubble collisions → turbulence) and the resulting GW spectrum, for a U(1)-charged waterfall field. |
| `dufaux2007_theory_numerics_gw_preheating_0707.0875.pdf` | Dufaux, Bergman, Felder, Kofman & Uzan, *Theory and Numerics of Gravitational Waves from Preheating after Inflation*, arXiv:0707.0875 | The other main early theory/numerics reference for computing Ω_GW(f) from classical field inhomogeneities at preheating — independent derivation to cross-check the CosmoLattice method against. |
| `gwabelian_cosmicstrings_preheating_1006.0217.pdf` | Dufaux, Figueroa & García-Bellido, *Gravitational Waves from Abelian Gauge Fields and Cosmic Strings at Preheating*, Phys. Rev. D 82, 083518 (2010), arXiv:1006.0217 | Directly on-topic: isolates what a U(1) gauge field (vs. scalar-only) adds to the preheating GW spectrum, including cosmic-string contributions — the closest precedent for this project's `lphi4` vs. `lphi4U1` comparison. |
| `tranberg2017_gw_nonabelian_tachyonic_1706.02365.pdf` | Tranberg, Tähtinen & Weir, *Gravitational waves from non-Abelian gauge fields at a tachyonic transition*, JCAP 04:012 (2018), arXiv:1706.02365 | SU(2)-Higgs (electroweak-like) tachyonic preheating GW spectrum — reports two spectral peaks tied to the Higgs and gauge-field masses. Directly relevant to the `lphi4SU2U1` side of this project. |
| `adshead2018_gw_gauge_preheating_1805.04550.pdf` | Adshead, Giblin & Weiner, *Gravitational waves from gauge preheating*, Phys. Rev. D 98, 043525 (2018), arXiv:1805.04550 | Modern reference on inflaton-gauge-field couplings (dilatonic and axial/Chern-Simons) driving efficient GW production during preheating — the comparison case for a gauge sector coupled through the inflaton rather than through the waterfall/Higgs field. |
| `gw_higgs_preheating_z2symmetry_2605.04670.pdf` | Zhou, Yu, Cheng & Zhang, *Gravitational Waves from Higgs Preheating after Inflaton Z2-Symmetry Breaking*, arXiv:2605.04670 (2026) | A recent CosmoLattice-based study (trilinear + quartic inflaton-Higgs couplings) showing current practice for extracting a GW spectrum with this exact code — a template for how to present and check the results this project will produce. |

## What's deliberately not here

General cosmology/inflation textbooks and broad GW-background reviews
aren't included — this folder is scoped to CosmoLattice itself and to
gauge fields at preheating specifically, not a general GW bibliography.
