"""Las figuras de la presentación con la frecuencia física de hoy en el eje x (en lugar de k) y h²Ω_GW hoy en el eje y.

Son las mismas curvas que las figuras del informe, con dos cambios de unidades que no tocan ninguna conclusión:
  f hoy = k · F_K,     F_K = 3.198×10⁷ Hz por unidad de k (data/produccion_gws.json, "f_por_unidad_de_k"; igual en los
                       tres modelos a 0.02 %, porque el traslado a hoy depende solo de la energía al final de la corrida)
  h²Ω_GW hoy = dΩ_GW/d ln k (final) · H2_FAC,  el mismo factor constante que usa produccion_gws.gws_hoy (g★ = 100)
(d ln k = d ln f, así que la densidad por intervalo logarítmico no cambia). Los cocientes entre modelos son idénticos.
La zona no convergida k > 4 pasa a ser f > 4·F_K ≈ 1.3×10⁸ Hz.

  figures/presentacion/espectro_final.png        <- figures/pendientes_espectro.png              (diapositiva 4)
  figures/presentacion/control_vs_U1.png         <- figures/comparacion_modelos/control_vs_U1.png (R2)
  figures/presentacion/ondas_gravitacionales.png <- .../lphi4SU2U1_vacioT_N64_kIR0.5_VV2/ondas_gravitacionales.png (R4)
  figures/presentacion/epocas_de_produccion.png  <- figures/produccion_gws/epocas_de_produccion.png (R5)

Uso: python3 code/figuras_presentacion.py   (después de code/produccion_gws.py y code/analisis_SU2U1.py)
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
import comparar_modelos  # noqa: E402
import pendientes_espectro  # noqa: E402
import produccion_gws as pg  # noqa: E402
from convergencia_N import C_AQUA, C_NARANJA, C_TINTA2, leer_espectros  # noqa: E402

DATA = FINAL / "data"
FIG = FINAL / "figures" / "presentacion"
F_K = json.load(open(DATA / "produccion_gws.json"))["hoy"]["control"]["f_por_unidad_de_k"]
H2_FAC = pg.H2_OMEGA_GAMMA * pg.GS0 ** (4 / 3) / 2 * pg.G_STAR ** (-1 / 3)
F_CONV = 4 * F_K                       # k = 4: borde de lo convergido
F_MAX = 30 * F_K
ETQ_F = "f hoy (Hz)"
ETQ_F7 = "f hoy (10⁷ Hz)"
ETQ_OM = "h² Ω_GW hoy"
NO_CONV = f"f > {F_CONV / 1e8:.1f}×10⁸ Hz: no convergido".replace(".", ",")


def positivo(y):
    return np.where(y > 0, y, np.nan)


def no_convergido(ax, texto=True):
    ax.axvspan(F_CONV, F_MAX, color="#eeeeea", zorder=0)
    if texto:
        ax.text(0.97, 0.05, NO_CONV, transform=ax.transAxes, ha="right", color=C_TINTA2, fontsize=8)


def espectro_final():
    fig, ax = plt.subplots(figsize=(7, 4.2))
    estilos = {1234: "-", 2345: "--", 3456: ":"}
    for m, (s, color, carpetas) in pendientes_espectro.CORRIDAS.items():
        for sem, c in carpetas.items():
            esp = leer_espectros(DATA / c / "spectra_energy_gws.txt")
            k, y = esp[-1, :, 0] / s, esp[-1, :, 1]
            sel = (k > 0) & (k <= 4 + 1e-3)
            ax.loglog(k[sel] * F_K, y[sel] * H2_FAC, color=color, ls=estilos[sem], lw=1.5, label=f"{m}, semilla {sem}")
    for (a, b), nombre in (((0.5, 1.5), "ajuste IR"), ((2.0, 4.0), "ajuste UV")):
        ax.axvspan(a * F_K, b * F_K, color="#eeeeea", zorder=0)
        ax.text(np.sqrt(a * b) * F_K, 0.03, nombre, color=C_TINTA2, fontsize=8, ha="center",
                transform=ax.get_xaxis_transform())
    ax.set_xlabel(ETQ_F)
    ax.set_ylabel(ETQ_OM)
    ax.set_title(f"Espectro de GWs hoy en la parte convergida (f ≤ {F_CONV / 1e8:.2f}×10⁸ Hz)".replace(".", ","))
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "espectro_final.png", dpi=150, bbox_inches="tight")


def control_vs_U1():
    R = [(et, comparar_modelos.cargar(c, s, cols), col) for et, c, s, cols, col in comparar_modelos.CORRIDAS]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for et, r, col in R:
        ax[0].semilogy(r["t"], np.maximum(r["frac"], 1e-16), color=col, label=et)
        ax[1].semilogy(r["t_gw"], positivo(r["rhoGW"]), color=col, label=et)
        ax[2].loglog(r["k"] * F_K, positivo(r["gw"]) * H2_FAC, color=col, label=et)
    ax[0].set_ylim(1e-14, 1)
    ax[0].set_ylabel("fracción de la energía en χ o en el gauge")
    ax[0].set_title("Energía transferida al campo hijo")
    ax[1].set_ylabel("ρ_GW / ρ")
    ax[1].set_title("Energía total en GWs")
    for a in ax[:2]:
        a.set_xlabel("τ (unidades del control)")
    no_convergido(ax[2])
    ax[2].set_xlabel(ETQ_F)
    ax[2].set_ylabel(ETQ_OM)
    ax[2].set_title("Espectro de GWs hoy")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "control_vs_U1.png", dpi=150, bbox_inches="tight")


def ondas_gravitacionales():
    res = json.load(open(DATA / "lphi4SU2U1_vacioT_N64_kIR0.5_VV2" / "analisis.json"))["gw"]
    R = {m: pg.cargar(p) for m, p in pg.MODELOS.items()}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for m, r in R.items():
        c = pg.MODELOS[m]["color"]
        ax[0].semilogy(r["t_gw"], positivo(r["rhoGW"]), color=c, lw=1.6, label=m)
        ax[1].loglog(r["k"] * F_K, positivo(r["om"][-1]) * H2_FAC, color=c, lw=1.6, label=m)
    for et, c, ls in (("SU2U1/control", C_NARANJA, "-"), ("SU2U1/U1", C_NARANJA, "--"), ("U1/control", C_AQUA, "-")):
        ax[2].plot(np.array(res[et]["k"]) * F_K / 1e7, res[et]["cociente"], color=c, ls=ls, lw=1.6, marker="o", ms=4,
                   label=et.replace("SU2U1", "SU(2)×U(1)").replace("U1", "U(1)"))
    ax[2].axhline(1, color=C_TINTA2, lw=0.8)
    no_convergido(ax[1])
    ax[0].set_xlabel("τ (unidades del control)")
    ax[0].set_ylabel("ρ_GW / ρ")
    ax[0].set_title("Energía total en GWs")
    ax[1].set_xlabel(ETQ_F)
    ax[1].set_ylabel(ETQ_OM)
    ax[1].set_title("Espectro de GWs hoy")
    ax[2].set_xlabel(ETQ_F7)
    ax[2].set_ylabel("cociente de espectros")
    ax[2].set_title(f"Cocientes en f = {0.5 * F_K / 1e7:.1f}–{2.5 * F_K / 1e7:.1f}×10⁷ Hz (una semilla)".replace(".", ","))
    ax[2].set_ylim(0, None)
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "ondas_gravitacionales.png", dpi=150, bbox_inches="tight")


def epocas_de_produccion():
    R = {m: pg.cargar(p) for m, p in pg.MODELOS.items()}
    k = R["control"]["k"]
    sel = (k > 0) & (k <= 6)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for m, r in R.items():
        c = pg.MODELOS[m]["color"]
        i_ep = int(np.argmin(abs(r["tau_esp"] - pg.EPOCA_W)))
        ax[0].plot(k[sel] * F_K / 1e7, r["om"][i_ep][sel] / r["om"][-1][sel], color=c, marker="o", ms=3.5, lw=1.5, label=m)
        acum = np.array([pg.integral_lnk(k, r["om"][i], 4) for i in range(len(r["tau_esp"]))]) / pg.integral_lnk(k, r["om"][-1], 4)
        ax[1].plot(r["tau_esp"], acum, color=c, lw=1.6, label=m)
    ax[0].axhline(1, color=C_TINTA2, lw=0.8)
    ax[0].axvspan(F_CONV / 1e7, 6 * F_K / 1e7, color="#eeeeea", zorder=0)
    ax[0].set_xlabel(ETQ_F7)
    ax[0].set_ylabel(f"Ω_GW(τ = {pg.EPOCA_W}) / Ω_GW(final)")
    ax[0].set_title(f"Qué parte de las ondas finales ya existía en τ = {pg.EPOCA_W}")
    ax[0].set_ylim(0, None)
    ax[1].axvspan(110, 130, color="#eeeeea", zorder=0)
    ax[1].text(120, 0.05, "disparo de\nlas W", ha="center", color=C_TINTA2, fontsize=8)
    ax[1].set_xlabel("τ (unidades del control)")
    ax[1].set_ylabel("Ω_GW acumulado / final")
    ax[1].set_title(f"Cuándo se producen las ondas (f ≤ {F_CONV / 1e8:.1f}×10⁸ Hz)".replace(".", ","))
    ax[1].set_xlim(40, 300)
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "epocas_de_produccion.png", dpi=150, bbox_inches="tight")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    espectro_final()
    control_vs_U1()
    ondas_gravitacionales()
    epocas_de_produccion()
    print(f"F_K = {F_K:.4e} Hz por unidad de k; H2_FAC = {H2_FAC:.4e}; f(k=4) = {F_CONV:.3e} Hz; f(k=2) = {2 * F_K:.3e} Hz")


if __name__ == "__main__":
    main()
