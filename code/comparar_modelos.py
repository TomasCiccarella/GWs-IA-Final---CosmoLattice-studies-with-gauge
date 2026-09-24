"""Compara el control lphi4 con los modelos gauge en la misma red (N = 64, kIR = 0.5 en unidades del control).

Corridas (las que existan y hayan terminado):
  control:                 data/convergencia_N/lphi4_N64_kIR0.5_VV2
  U(1), CosmoLattice:      data/lphi4U1_N64_kIR0.5_VV2         (A = 0, solo E longitudinal de Gauss)
  U(1), vacío transversal: data/lphi4U1_vacioT_N64_kIR0.5_VV2  (parche code/parches/u1_vacio_transversal.patch)

Todo en unidades del control: k = k_gauge/√2, τ = √2 t_gauge. dΩ_GW/d ln k y ρ_GW/ρ son adimensionales.
Figuras en figures/comparacion_modelos/, números en data/comparacion_modelos.json.

Uso: python3 code/comparar_modelos.py
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from convergencia_N import (C_AQUA, C_NARANJA, C_TINTA2, leer_espectros, leer_promedios,  # noqa: E402
                            primer_cruce, terminada)

DATA = FINAL / "data"
FIG = FINAL / "figures" / "comparacion_modelos"
S2 = np.sqrt(2)
C_AZUL = "#2a78d6"

# (etiqueta, carpeta, escala de k y de t hacia unidades del control, columnas de energía del "hijo", color)
CORRIDAS = [
    ("control λφ⁴ + χ", "convergencia_N/lphi4_N64_kIR0.5_VV2", 1.0, ["E^kin_scal1", "E^grad_scal1"], C_AZUL),
    ("U(1), condición de CosmoLattice", "lphi4U1_N64_kIR0.5_VV2", S2, ["E^kin_U10", "E^grad_U10"], C_NARANJA),
    ("U(1), vacío transversal", "lphi4U1_vacioT_N64_kIR0.5_VV2", S2, ["E^kin_U10", "E^grad_U10"], C_AQUA),
]


def cargar(carpeta, s, cols):
    D = DATA / carpeta
    en = leer_promedios(D / "average_energies.txt")
    gw = leer_promedios(D / "average_energies_gws.txt")
    esp = leer_espectros(D / "spectra_energy_gws.txt")
    t_esp = np.atleast_1d(np.loadtxt(D / "average_spectra_times.txt"))[: len(esp)] * s
    frac = sum(en[c] for c in cols) / en["E_tot"]
    return dict(t=en["t"] * s, frac=frac, t_gw=gw["t"] * s, rhoGW=gw["rhoGW_over_rho"],
                t_esp=t_esp, k=esp[-1, :, 0] / s, gw=esp[-1, :, 1], esp=esp[:, :, 1])


def main():
    R = {et: (cargar(c, s, cols), col) for et, c, s, cols, col in CORRIDAS if terminada(DATA / c)}
    FIG.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for et, (r, col) in R.items():
        ax[0].semilogy(r["t"], np.maximum(r["frac"], 1e-16), color=col, label=et)
        ax[1].semilogy(r["t_gw"], np.where(r["rhoGW"] > 0, r["rhoGW"], np.nan), color=col, label=et)
        ax[2].loglog(r["k"], np.where(r["gw"] > 0, r["gw"], np.nan), color=col, label=et)
    ax[0].set_ylim(1e-14, 1)
    ax[0].set_ylabel("fracción de la energía en χ o en el gauge")
    ax[0].set_title("Energía transferida al campo hijo")
    ax[1].set_ylabel("ρ_GW / ρ")
    ax[1].set_title("Energía total en GWs")
    for a in ax[:2]:
        a.set_xlabel("τ (unidades del control)")
    ax[2].axvspan(4, 30, color="#eeeeea", zorder=0)
    ax[2].text(0.97, 0.05, "k > 4: no convergido", transform=ax[2].transAxes, ha="right",
               color=C_TINTA2, fontsize=8)
    ax[2].set_xlabel("k (unidades del control)")
    ax[2].set_ylabel("dΩ_GW / d ln k (final)")
    ax[2].set_title("Espectro de GWs al final (τ = 300)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "control_vs_U1.png", dpi=130, bbox_inches="tight")

    res = {}
    for et, (r, _) in R.items():
        jp = int(np.argmax(np.where(r["k"] <= 4, r["gw"], 0)))
        res[et] = dict(
            tau_10pc=primer_cruce(r["t"], r["frac"], 0.1), tau_1pc=primer_cruce(r["t"], r["frac"], 0.01),
            frac_max=float(r["frac"].max()), rhoGW_fin=float(r["rhoGW"][-1]),
            k_pico_k_menor_4=float(r["k"][jp]), gw_pico_k_menor_4=float(r["gw"][jp]))
    json.dump(res, open(DATA / "comparacion_modelos.json", "w"), indent=2, ensure_ascii=False)
    for et, d in res.items():
        print(f"{et:34s} τ(1 %) = {d['tau_1pc']:6.1f}  τ(10 %) = {d['tau_10pc']:6.1f}  frac máx = {d['frac_max']:.3f}  "
              f"ρGW/ρ fin = {d['rhoGW_fin']:.2e}  pico (k≤4): {d['gw_pico_k_menor_4']:.2e} en k = {d['k_pico_k_menor_4']:.2f}")


if __name__ == "__main__":
    main()
