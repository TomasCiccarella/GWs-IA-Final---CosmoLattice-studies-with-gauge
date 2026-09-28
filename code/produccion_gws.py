"""Cómo se producen las GWs en los tres modelos: cuatro análisis sobre las corridas que ya existen.

Corridas (N = 64, kIR = 0.5 en unidades del control, VV2, semilla 1234; unidades del control k = k_gauge/√2,
τ = √2 t_gauge):
  control:     data/convergencia_N/lphi4_N64_kIR0.5_VV2
  U(1):        data/lphi4U1_vacioT_N64_kIR0.5_VV2
  SU(2)×U(1):  data/lphi4SU2U1_vacioT_N64_kIR0.5_VV2

1. Épocas: el espectro dΩ_GW/d ln k cada 5 τ. Como el fondo se comporta como radiación (λφ⁴) y las GWs libres
   también, dΩ_GW/d ln k deja de cambiar cuando ya no hay fuente; así que Ω(τ)/Ω(final) por bin dice qué parte de
   las ondas finales ya estaba producida en τ. En SU(2)×U(1) las W se disparan en τ ≈ 110-130 (bitácora 3.7):
   comparar lo producido antes y después separa lo que aportan la Z (y el fotón) de lo que aportan las W.
2. Ritmo: en la fase lineal, si el campo crece como e^{μτ} su energía crece como e^{2μτ} y la fuente de las GWs
   (cuadrática en el campo) como e^{4μτ}, así que se espera d ln ρ_GW/dτ ≈ 2 d ln ρ_campo/dτ.
3. Eficiencia: Ω_GW (integrado en la parte convergida, k ≤ 2.5) dividido por el cuadrado de la fracción de la
   energía en el campo fuente durante la producción (τ = 70-150), que es la escala de la estimación estándar
   Ω_GW ~ (ρ_fuente/ρ)² (a H / k)²; el factor (aH/k)² es el mismo en los tres modelos (mismo pico k = 2).
4. Hoy: f y h²Ω_GW hoy con la misma cuenta que el piloto (code/analisis_corridas.py: radiación y conservación de
   la entropía desde el final de la corrida, g★ = 100).

Uso: python3 code/produccion_gws.py   (figuras en figures/produccion_gws/, números en data/produccion_gws.json)
"""
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from convergencia_N import C_AQUA, C_NARANJA, C_TINTA2, leer_espectros, leer_promedios  # noqa: E402

DATA = FINAL / "data"
FIG = FINAL / "figures" / "produccion_gws"
S2 = np.sqrt(2)
C_AZUL = "#2a78d6"
LAM = 9e-14
G_STAR = 100
K_CONV = 2.5            # parte convergida del espectro (bitácora 3.3)
EPOCA_W = 100           # τ antes del disparo de las W en SU(2)×U(1) (1 % en τ = 117)
VENTANA_LINEAL = (35, 62)
VENTANA_PRODUCCION = (70, 150)
TIEMPOS_FIG = [70, 85, 100, 115, 130, 150, 200, 295]

MODELOS = {
    "control": dict(carpeta="convergencia_N/lphi4_N64_kIR0.5_VV2", s=1.0, color=C_AZUL,
                    fuente=lambda e: e["E^kin_scal1"] + e["E^grad_scal1"], crece=lambda e: e["E^grad_scal1"],
                    fstar="initial_amplitudes"),
    "U(1)": dict(carpeta="lphi4U1_vacioT_N64_kIR0.5_VV2", s=S2, color=C_AQUA,
                 fuente=lambda e: e["E^kin_U10"] + e["E^grad_U10"], crece=lambda e: e["E^grad_U10"],
                 fstar="cmplx_field_initial_norm"),
    "SU(2)×U(1)": dict(carpeta="lphi4SU2U1_vacioT_N64_kIR0.5_VV2", s=S2, color=C_NARANJA,
                       fuente=lambda e: e["E^kin_U10"] + e["E^grad_U10"] + e["E^kin_SU2"] + e["E^grad_SU2"],
                       crece=lambda e: 2 * e["E^grad_U10"], fstar="SU2Doublet_initial_norm"),
}

# Constantes del traslado a hoy (las mismas que code/analisis_corridas.py)
HBAR_GEV_S = 6.582119569e-25
T0_GEV = 2.7255 * 8.617333262e-14
GS0 = 3.931
H2_OMEGA_GAMMA = 2.473e-5


def parametro(D, nombre):
    txt = next(D.glob("*.in")).read_text()
    return float(re.search(rf"^\s*{nombre}\s*=\s*(\S+)", txt, re.M).group(1))


def cargar(p):
    D = DATA / p["carpeta"]
    en = leer_promedios(D / "average_energies.txt")
    gw = leer_promedios(D / "average_energies_gws.txt")
    fondo = leer_promedios(D / "average_scale_factor.txt")
    esp = leer_espectros(D / "spectra_energy_gws.txt")
    t_esp = np.atleast_1d(np.loadtxt(D / "average_spectra_times.txt"))[: len(esp)]
    fstar = parametro(D, p["fstar"])
    return dict(t=en["t"] * p["s"], en=en, T=en["E_tot"], t_gw=gw["t"] * p["s"], rhoGW=gw["rhoGW_over_rho"],
                tau_esp=t_esp * p["s"], a_esp=np.interp(t_esp, fondo["t"], fondo["a"]),
                k_prog=esp[0, :, 0], k=esp[0, :, 0] / p["s"], om=esp[:, :, 1],
                fstar=fstar, omstar=np.sqrt(LAM) * fstar)


def pendiente(t, y, v):
    m = (t >= v[0]) & (t <= v[1]) & (y > 0)
    return float(np.polyfit(t[m], np.log(y[m]), 1)[0])


def integral_lnk(k, om, kmax):
    """∫ dΩ/d ln k d ln k sobre los bins con k ≤ kmax (bins lineales de ancho Δk: d ln k ≈ Δk/k)."""
    dk = k[1] - k[0]
    m = (k > 0) & (k <= kmax + 1e-3)   # tolerancia: en los modelos gauge k = n·0.5 sale con error de 1e-6
    return float(np.sum(om[m] * dk / k[m]))


def gws_hoy(r, g_star=G_STAR):
    rho_fin = r["T"][-1] * r["fstar"] ** 2 * r["omstar"] ** 2
    T_fin = (30 * rho_fin / (np.pi ** 2 * g_star)) ** 0.25
    a_fin_sobre_a0 = (GS0 / g_star) ** (1 / 3) * T0_GEV / T_fin
    f0 = r["k_prog"] * r["omstar"] / r["a_esp"][-1] * a_fin_sobre_a0 / (2 * np.pi * HBAR_GEV_S)
    h2 = H2_OMEGA_GAMMA * GS0 ** (4 / 3) / 2 * g_star ** (-1 / 3) * r["om"][-1]
    return f0, h2


def main():
    R = {m: cargar(p) for m, p in MODELOS.items()}
    k = R["control"]["k"]
    for r in R.values():
        assert np.allclose(r["k"], k)
    conv = (k > 0) & (k <= K_CONV + 1e-3)
    res = {"epocas": {}, "ritmo": {}, "eficiencia": {}, "hoy": {}}

    # 1. Épocas
    for m, r in R.items():
        i_ep = int(np.argmin(abs(r["tau_esp"] - EPOCA_W)))
        i_120 = int(np.argmin(abs(r["tau_esp"] - 120)))
        fin = r["om"][-1]
        res["epocas"][m] = dict(
            tau_epoca=float(r["tau_esp"][i_ep]),
            fraccion_producida_antes_por_bin=dict(zip(np.round(k[conv], 2).tolist(),
                                                      np.round(r["om"][i_ep][conv] / fin[conv], 3).tolist())),
            fraccion_antes_k_hasta_4=integral_lnk(k, r["om"][i_ep], 4) / integral_lnk(k, fin, 4),
            fraccion_antes_120_k_hasta_4=integral_lnk(k, r["om"][i_120], 4) / integral_lnk(k, fin, 4),
            Omega_convergido_antes=integral_lnk(k, r["om"][i_ep], K_CONV),
            tau_mitad_k_hasta_4=float(np.interp(0.5, np.maximum.accumulate(
                [integral_lnk(k, r["om"][i], 4) / integral_lnk(k, fin, 4) for i in range(len(r["tau_esp"]))]),
                r["tau_esp"])),
            rhoGW_tau75=float(np.interp(75, r["t_gw"], r["rhoGW"])),
            rhoGW_tau60=float(np.interp(60, r["t_gw"], r["rhoGW"])),
            Omega_convergido_final=integral_lnk(k, fin, K_CONV))
    # Lo que produjo SU(2)×U(1) antes de τ = 100 comparado con el control y U(1) en el mismo instante
    for m in ("control", "U(1)"):
        res["epocas"][f"SU2U1/{m} antes de τ={EPOCA_W}"] = (res["epocas"]["SU(2)×U(1)"]["Omega_convergido_antes"] /
                                                          res["epocas"][m]["Omega_convergido_antes"])
        res["epocas"][f"SU2U1/{m} final"] = (res["epocas"]["SU(2)×U(1)"]["Omega_convergido_final"] /
                                             res["epocas"][m]["Omega_convergido_final"])

    # 2. Ritmo de crecimiento en la fase lineal
    for m, r in R.items():
        s_campo = pendiente(r["t"], MODELOS[m]["crece"](r["en"]) / r["T"], VENTANA_LINEAL)
        s_gw = pendiente(r["t_gw"], r["rhoGW"], VENTANA_LINEAL)
        res["ritmo"][m] = dict(dlnrho_campo=s_campo, dlnrho_gw=s_gw, cociente=s_gw / s_campo,
                               mu_campo=s_campo / 2, mu_gw_equivalente=s_gw / 4)

    # 3. Eficiencia
    for m, r in R.items():
        f = MODELOS[m]["fuente"](r["en"]) / r["T"]
        v = (r["t"] >= VENTANA_PRODUCCION[0]) & (r["t"] <= VENTANA_PRODUCCION[1])
        f_media = float(f[v].mean())
        om_c = res["epocas"][m]["Omega_convergido_final"]
        res["eficiencia"][m] = dict(fraccion_media_produccion=f_media, fraccion_final=float(f[-1]),
                                    Omega_convergido=om_c, eficiencia=om_c / f_media ** 2)
    for m in ("U(1)", "SU(2)×U(1)"):
        res["eficiencia"][f"{m}/control"] = dict(
            Omega=res["eficiencia"][m]["Omega_convergido"] / res["eficiencia"]["control"]["Omega_convergido"],
            fraccion_al_cuadrado=(res["eficiencia"][m]["fraccion_media_produccion"] /
                                  res["eficiencia"]["control"]["fraccion_media_produccion"]) ** 2,
            eficiencia=res["eficiencia"][m]["eficiencia"] / res["eficiencia"]["control"]["eficiencia"])

    # 4. Hoy
    HOY = {}
    for m, r in R.items():
        f0, h2 = gws_hoy(r)
        HOY[m] = (f0, h2)
        jp = int(np.argmax(np.where(k <= 4, h2, 0)))
        jg = int(np.argmax(h2))
        res["hoy"][m] = dict(f_pico_k_menor_4=float(f0[jp]), h2Omega_pico_k_menor_4=float(h2[jp]),
                             f_pico_global=float(f0[jg]), h2Omega_pico_global=float(h2[jg]),
                             f_por_unidad_de_k=float(f0[1] / k[1]))
    json.dump(res, open(DATA / "produccion_gws.json", "w"), indent=2, ensure_ascii=False)

    # ---------------- Figuras ----------------
    FIG.mkdir(parents=True, exist_ok=True)
    norm = Normalize(min(TIEMPOS_FIG), max(TIEMPOS_FIG))
    cmap = plt.get_cmap("viridis")

    # 1a. Espectros en el tiempo
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for j, (m, r) in enumerate(R.items()):
        for tau in TIEMPOS_FIG:
            i = int(np.argmin(abs(r["tau_esp"] - tau)))
            y = r["om"][i]
            ax[j].loglog(k[1:], np.where(y[1:] > 0, y[1:], np.nan), color=cmap(norm(tau)), lw=1.3,
                         label=f"τ = {r['tau_esp'][i]:.0f}")
        ax[j].axvspan(4, 30, color="#eeeeea", zorder=0)
        ax[j].set_title(m)
        ax[j].set_xlabel("k (unidades del control)")
        ax[j].set_ylim(1e-10, 3e-5)
    ax[0].set_ylabel("dΩ_GW / d ln k")
    ax[2].legend(fontsize=7, loc="lower left", title="espectro en", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "espectros_en_el_tiempo.png", dpi=130, bbox_inches="tight")

    # 1b. Fracción de las GWs finales producida antes de τ = 100, por bin
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    kk = k[(k > 0) & (k <= 6)]
    for m, r in R.items():
        i_ep = int(np.argmin(abs(r["tau_esp"] - EPOCA_W)))
        sel = (k > 0) & (k <= 6)
        ax[0].plot(kk, r["om"][i_ep][sel] / r["om"][-1][sel], color=MODELOS[m]["color"], marker="o", ms=3.5, lw=1.5,
                   label=m)
        t_all = r["tau_esp"]
        acum = np.array([integral_lnk(k, r["om"][i], 4) for i in range(len(t_all))]) / integral_lnk(k, r["om"][-1], 4)
        ax[1].plot(t_all, acum, color=MODELOS[m]["color"], lw=1.6, label=m)
    ax[0].axhline(1, color=C_TINTA2, lw=0.8)
    ax[0].axvspan(4, 6, color="#eeeeea", zorder=0)
    ax[0].set_xlabel("k (unidades del control)")
    ax[0].set_ylabel(f"Ω_GW(τ = {EPOCA_W}) / Ω_GW(final)")
    ax[0].set_title(f"Qué parte de las ondas finales ya existía en τ = {EPOCA_W}")
    ax[0].set_ylim(0, None)
    ax[1].axvspan(110, 130, color="#eeeeea", zorder=0)
    ax[1].text(120, 0.05, "disparo de\nlas W", ha="center", color=C_TINTA2, fontsize=8)
    ax[1].set_xlabel("τ (unidades del control)")
    ax[1].set_ylabel("Ω_GW(k ≤ 4) acumulado / final")
    ax[1].set_title("Cuándo se producen las ondas (k ≤ 4)")
    ax[1].set_xlim(40, 300)
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "epocas_de_produccion.png", dpi=130, bbox_inches="tight")

    # 2. Ritmo
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for m, r in R.items():
        c = MODELOS[m]["color"]
        campo = MODELOS[m]["crece"](r["en"]) / r["T"]
        ax.semilogy(r["t"], np.where(campo > 0, campo, np.nan) ** 2, color=c, lw=1.0, alpha=0.6)
        ax.semilogy(r["t_gw"], np.where(r["rhoGW"] > 0, r["rhoGW"], np.nan), color=c, lw=2.0,
                    label=f"{m}: ρ_GW/ρ (cociente de pendientes {res['ritmo'][m]['cociente']:.2f})")
    ax.axvspan(*VENTANA_LINEAL, color="#eeeeea", zorder=0)
    ax.plot([], [], color=C_TINTA2, lw=1.0, alpha=0.6, label="(ρ_campo/ρ)², misma escala de color")
    ax.set_xlim(0, 120)
    ax.set_ylim(1e-30, 1e-3)
    ax.set_xlabel("τ (unidades del control)")
    ax.set_title("La energía en GWs crece como el cuadrado de la del campo")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "ritmo_gws_vs_campo.png", dpi=130, bbox_inches="tight")

    # 4. Hoy
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for m, (f0, h2) in HOY.items():
        ok = (h2 > 0) & (k > 0)
        ax.loglog(f0[ok], h2[ok], color=MODELOS[m]["color"], lw=1.8, label=m)
        conv_ok = ok & (k > 4)
    f_k4 = HOY["control"][0][np.argmin(abs(k - 4))]
    ax.axvspan(f_k4, HOY["control"][0][-1] * 1.2, color="#eeeeea", zorder=0)
    ax.text(0.97, 0.05, "k > 4: no convergido", transform=ax.transAxes, ha="right", color=C_TINTA2, fontsize=8)
    ax.set_xlabel("f hoy (Hz)")
    ax.set_ylabel("h² Ω_GW hoy")
    ax.set_title(f"Espectro de GWs hoy (g★ = {G_STAR})")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "espectro_hoy.png", dpi=130, bbox_inches="tight")

    # Resumen
    print(json.dumps({x: res[x] for x in ("ritmo", "eficiencia", "hoy")}, indent=1, ensure_ascii=False))
    for m in R:
        e = res["epocas"][m]
        print(m, "fracción producida antes de τ=100 (k≤4):", round(e["fraccion_antes_k_hasta_4"], 3),
              "antes de 120:", round(e["fraccion_antes_120_k_hasta_4"], 3), "τ mitad:", round(e["tau_mitad_k_hasta_4"], 1),
              "ρGW τ60/75:", f"{e['rhoGW_tau60']:.1e} {e['rhoGW_tau75']:.1e}", "por bin:", e["fraccion_producida_antes_por_bin"])
    print({x: y for x, y in res["epocas"].items() if x.startswith("SU2U1/")})


if __name__ == "__main__":
    main()
