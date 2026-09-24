"""Crecimiento por modo del campo gauge en lphi4U1 contra Floquet y contra χ del control.

Mide μ_k = d(½ ln Δ)/dτ en ventanas de la fase lineal, para el campo magnético (transversal) U(1) de
data/lphi4U1_N64_kIR0.5_VV2 y para χ del control con la misma red
(data/convergencia_N/lphi4_N64_kIR0.5_VV2). Todo en unidades del control:
k_lphi4 = k_gauge/√2 y τ = √2 t_gauge (ver analisis_parametros.traducir).
Floquet: ecuación de Lamé con amplitud conforme A = 1.138 (medida en el piloto), para q = 120
(el emparejado) y q = 60 y 240 como diagnóstico de un factor 2.

Uso: python3 code/crecimiento_por_modo_U1.py [carpeta en data/, por defecto lphi4U1_N64_kIR0.5_VV2]
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from convergencia_N import C_AQUA, C_NARANJA, C_TINTA2, leer_espectros, leer_promedios  # noqa: E402
from analisis_parametros import floquet_mu  # noqa: E402

U1 = FINAL / "data" / (sys.argv[1] if len(sys.argv) > 1 else "lphi4U1_N64_kIR0.5_VV2")
CONTROL = FINAL / "data" / "convergencia_N" / "lphi4_N64_kIR0.5_VV2"
FIG = FINAL / "figures" / U1.name
S2 = np.sqrt(2)
A_CONF = 1.138
VENTANAS = [(10, 30), (30, 50), (50, 70)]   # en τ del control
C_AZUL = "#2a78d6"


def cargar(D, archivo, escala_k, escala_t, col=1):
    """col = 1: primer espectro del archivo (Δ del campo para χ; campo ELÉCTRICO para U(1));
    col = 2: segundo (Δ del momento para χ; campo MAGNÉTICO para U(1), ver u1measurer.h)."""
    esp = leer_espectros(D / archivo)
    t = np.atleast_1d(np.loadtxt(D / "average_spectra_times.txt"))[: len(esp)] * escala_t
    fondo = leer_promedios(D / "average_scale_factor.txt")
    a = np.interp(t / escala_t, fondo["t"], fondo["a"])
    return esp[0, :, 0] * escala_k, t, a, esp[:, :, col]


def mu_por_modo(t, a, delta, ventana, potencia_a):
    """Pendiente de ½ ln(a^p Δ) en la ventana (p = 2 para un escalar, 0 para el gauge, que es conforme)."""
    m = (t >= ventana[0]) & (t <= ventana[1])
    lnX = 0.5 * np.log(a[:, None] ** potencia_a * np.maximum(delta, 1e-300))
    return np.polyfit(t[m], lnX[m], 1)[0]


def main():
    # Campo magnético: es puramente transversal (B = rot A), así que mide los modos que resuenan. El eléctrico
    # incluye la parte longitudinal fijada por Gauss, que no resuena y tapa el crecimiento.
    kU, tU, aU, dU = cargar(U1, "spectra_norm_U1_0.txt", 1 / S2, S2, col=2)
    kC, tC, aC, dC = cargar(CONTROL, "spectra_scalar_1.txt", 1.0, 1.0)
    Kf = np.linspace(0.02, 4, 120)
    floq = {q: np.array([A_CONF * floquet_mu(K / A_CONF, q, pasos=1500) for K in Kf]) for q in (60, 120, 240)}

    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, len(VENTANAS), figsize=(13, 3.8), sharey=True)
    res = {}
    for j, v in enumerate(VENTANAS):
        if tU[-1] < v[1]:
            continue
        muU = mu_por_modo(tU, aU, dU, v, 0)
        muC = mu_por_modo(tC, aC, dC, v, 2)
        for q, ls in ((120, "-"), (60, ":"), (240, "--")):
            ax[j].plot(Kf, floq[q], color=C_AZUL, ls=ls, lw=1.2 if q == 120 else 0.9, label=f"Floquet q = {q}")
        ax[j].plot(kC, muC, "o", color=C_AQUA, ms=4, mfc="white", mew=1.3, label="χ del control (misma red)")
        ax[j].plot(kU, muU, "s", color=C_NARANJA, ms=4, label="campo U(1) (espectro magnético)")
        ax[j].axhline(0, color=C_TINTA2, lw=0.6)
        ax[j].set_xlim(0, 4)
        ax[j].set_title(f"τ ∈ [{v[0]}, {v[1]}]")
        ax[j].set_xlabel("k (unidades del control)")
        enb = (kU > 0.4) & (kU < 2.2)
        res[f"{v[0]}-{v[1]}"] = dict(
            U1_mu_max=float(np.max(muU[enb])), U1_k_max=float(kU[enb][np.argmax(muU[enb])]),
            chi_mu_max=float(np.max(muC[(kC > 0.4) & (kC < 2.2)])),
            U1=dict(zip(np.round(kU[kU < 4], 3).tolist(), np.round(muU[kU < 4], 4).tolist())))
    ax[0].set_ylabel("μ_k (por unidad de τ)")
    ax[0].set_ylim(-0.1, 0.35)
    ax[0].legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "mu_por_modo_vs_control.png", dpi=130, bbox_inches="tight")
    res["floquet_q120_max"] = dict(mu=float(floq[120].max()), k=float(Kf[np.argmax(floq[120])]))
    json.dump(res, open(U1 / "crecimiento_por_modo.json", "w"), indent=2)
    for v, r in res.items():
        if "U1_mu_max" in r:
            print(f"τ {v}: U(1) μ_max = {r['U1_mu_max']:.3f} en k = {r['U1_k_max']:.2f};  χ control μ_max = {r['chi_mu_max']:.3f}")
    print("Floquet q=120:", res["floquet_q120_max"])


if __name__ == "__main__":
    main()
