"""Convergencia en N del control lphi4 (VV2): ¿cuánto cambian los espectros si bajamos la resolución?

Compara las corridas de data/convergencia_N/ (N = 32, 48, 64, 96 con kIR = 0.25, y N = 64 con
kIR = 0.5) contra la referencia data/lphi4_control_N128_VV2. Todas tienen el mismo .in salvo N
(y kIR en la última). Guarda figuras en figures/convergencia_N/ y los números en
data/convergencia_N/resumen.json.

Uso: python3 code/convergencia_N.py
"""
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

FINAL = Path(__file__).resolve().parent.parent
DATA = FINAL / "data"
FIG = FINAL / "figures" / "convergencia_N"

REFERENCIA = "lphi4_control_N128_VV2"
CORRIDAS = [  # (carpeta, etiqueta)
    ("convergencia_N/lphi4_N32_VV2", "N = 32"),
    ("convergencia_N/lphi4_N48_VV2", "N = 48"),
    ("convergencia_N/lphi4_N64_VV2", "N = 64"),
    ("convergencia_N/lphi4_N96_VV2", "N = 96"),
    (REFERENCIA, "N = 128 (referencia)"),
    ("convergencia_N/lphi4_N64_kIR0.5_VV2", "N = 64, kIR = 0.5"),
    ("convergencia_N/lphi4_N128_kIR0.5_VV2", "N = 128, kIR = 0.5"),
]

# Paleta del notebook de análisis: rampa azul para N (claro = N chico), naranja/aqua para kIR = 0.5 (N = 64/128)
C_NARANJA, C_AQUA, C_TINTA, C_TINTA2 = "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
RAMPA = LinearSegmentedColormap.from_list("rampa_azul", ["#a9c6ea", "#2a78d6", "#0b2a55"])
plt.rcParams.update({
    "figure.dpi": 110, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
    "axes.edgecolor": "#b9b8b2", "axes.labelcolor": C_TINTA, "xtick.color": C_TINTA2,
    "ytick.color": C_TINTA2, "lines.linewidth": 1.5, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "font.size": 10,
    "axes.titlesize": 10, "axes.titlelocation": "left",
})

# Constantes para llevar el espectro a hoy (mismas que code/analisis_corridas.py)
HBAR_GEV_S, KB_GEV_K, GS0, H2_OMEGA_GAMMA, G_STAR = 6.582119569e-25, 8.617333262e-14, 3.931, 2.473e-5, 100
T0_GEV = 2.7255 * KB_GEV_K


def leer_promedios(ruta):
    with open(ruta) as fh:
        cabecera = fh.readline().lstrip("#").split()
    datos = np.loadtxt(ruta, ndmin=2)
    return {nombre: datos[:, i] for i, nombre in enumerate(cabecera)}


def leer_espectros(ruta):
    bloques, actual = [], []
    for linea in open(ruta):
        linea = linea.strip()
        if not linea:
            if actual:
                bloques.append(actual)
                actual = []
            continue
        if ":" in linea or linea.startswith("#"):
            continue
        actual.append([float(x) for x in linea.split()])
    if actual:
        bloques.append(actual)
    return np.array(bloques)


def leer_parametros(ruta_in):
    pars = {}
    for linea in open(ruta_in):
        linea = linea.split("#")[0].strip()
        if "=" in linea:
            clave, valor = (s.strip() for s in linea.split("=", 1))
            pars[clave] = valor
    return pars


def primer_cruce(t, y, umbral):
    i = np.flatnonzero(y >= umbral)
    if len(i) == 0:
        return np.nan
    i = i[0]
    return t[0] if i == 0 else float(np.interp(np.log(umbral), np.log([y[i - 1], y[i]]), [t[i - 1], t[i]]))


def leer_tiempo(ruta):
    """Tiempo de pared (s) y memoria máxima (MB) de un tiempo.log de /usr/bin/time -v."""
    txt = open(ruta).read()
    h = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", txt).group(1).split(":")
    seg = sum(float(x) * 60**i for i, x in enumerate(reversed(h)))
    mem = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1)) / 1024
    return seg, mem


def terminada(D):
    return (D / "tiempo.log").exists() and "Exit status: 0" in (D / "tiempo.log").read_text()


def cargar(carpeta):
    D = DATA / carpeta
    P = leer_parametros(sorted(D.glob("*.in"))[0])
    fondo = leer_promedios(D / "average_scale_factor.txt")
    en = leer_promedios(D / "average_energies.txt")
    gw = leer_promedios(D / "average_energies_gws.txt")
    fr = leer_promedios(D / "average_energy_conservation.txt")
    esp_chi = leer_espectros(D / "spectra_scalar_1.txt")
    esp_phi = leer_espectros(D / "spectra_scalar_0.txt")
    esp_gw = leer_espectros(D / "spectra_energy_gws.txt")
    t_esp = np.atleast_1d(np.loadtxt(D / "average_spectra_times.txt"))[: len(esp_gw)]
    a_fin = float(np.interp(t_esp[-1], fondo["t"], fondo["a"]))
    lam, fStar = float(P["lambda"]), float(P["initial_amplitudes"].split()[0])
    omegaStar = np.sqrt(lam) * fStar

    # Espectro de GWs hoy (radiación hasta hoy, g★ = 100), igual que el notebook
    k = esp_gw[-1, :, 0]
    rho_fin = float(en["E_tot"][-1]) * fStar**2 * omegaStar**2
    T_fin = (30 * rho_fin / (np.pi**2 * G_STAR)) ** 0.25
    f0 = k * omegaStar / a_fin * (GS0 / G_STAR) ** (1 / 3) * T0_GEV / T_fin / (2 * np.pi * HBAR_GEV_S)
    h2O = H2_OMEGA_GAMMA * GS0 ** (4 / 3) / 2 * G_STAR ** (-1 / 3) * esp_gw[-1, :, 1]

    # Fracción de la energía que está en χ (cinética + gradiente): mide cuándo termina la resonancia
    f_chi = (en["E^kin_scal1"] + en["E^grad_scal1"]) / en["E_tot"]
    jp = int(np.argmax(esp_gw[-1, :, 1]))
    seg, mem = leer_tiempo(D / "tiempo.log")
    return dict(
        N=int(P["N"]), kIR=float(P["kIR"]), t=en["t"], f_chi=f_chi, t_gw=gw["t"], rhoGW=gw["rhoGW_over_rho"],
        k_chi=esp_chi[-1, :, 0], chi=esp_chi[-1, :, 1], k_phi=esp_phi[-1, :, 0], phi=esp_phi[-1, :, 1],
        k=k, gw=esp_gw[-1, :, 1], f0=f0, h2O=h2O,
        tau_res=primer_cruce(en["t"], f_chi, 0.1), rhoGW_fin=float(gw["rhoGW_over_rho"][-1]),
        k_pico=float(k[jp]), gw_pico=float(esp_gw[-1, jp, 1]), f_pico=float(f0[jp]), h2O_pico=float(h2O[jp]),
        friedmann_max=float(np.max(np.abs(fr["rel_diff_friedmann"]))), segundos=seg, memoria_MB=mem,
        a_fin=a_fin,
    )


def en_comun(kx, yx, kref, yref):
    """Cociente y/y_ref sobre los bins de k que comparten (los bins son múltiplos de kIR)."""
    comunes, ix, ir = np.intersect1d(np.round(kx, 8), np.round(kref, 8), return_indices=True)
    ok = (yref[ir] > 0) & (yx[ix] > 0)
    return comunes[ok], yx[ix][ok] / yref[ir][ok]


def main():
    R = {et: cargar(c) for c, et in CORRIDAS if terminada(DATA / c)}
    ref = R["N = 128 (referencia)"]
    estilo = {}
    Ns = [r["N"] for et, r in R.items() if r["kIR"] == 0.25]
    for et, r in R.items():
        if r["kIR"] == 0.25:
            x = (np.log(r["N"]) - np.log(min(Ns))) / max(np.log(max(Ns)) - np.log(min(Ns)), 1e-9)
            estilo[et] = dict(color=RAMPA(x), ls="-", lw=2.2 if r is ref else 1.4)
        else:
            estilo[et] = dict(color=C_NARANJA if r["N"] == 64 else C_AQUA, ls="--", lw=1.4)
    FIG.mkdir(parents=True, exist_ok=True)

    # 1) Espectro de GWs al final y cociente contra N = 128
    fig, ax = plt.subplots(2, 2, figsize=(12, 6.5), sharex="col", gridspec_kw=dict(height_ratios=[2.2, 1]))
    for et, r in R.items():
        ax[0, 0].loglog(r["k"], np.where(r["gw"] > 0, r["gw"], np.nan), label=et, **estilo[et])
        ax[0, 1].loglog(r["k_chi"], r["chi"], label=et, **estilo[et])
        if r is not ref:
            ax[1, 0].semilogx(*en_comun(r["k"], r["gw"], ref["k"], ref["gw"]), **estilo[et])
            ax[1, 1].semilogx(*en_comun(r["k_chi"], r["chi"], ref["k_chi"], ref["chi"]), **estilo[et])
    ax[0, 0].set_ylim(ref["gw_pico"] * 1e-6, max(r["gw_pico"] for r in R.values()) * 3)
    ax[0, 0].set_ylabel("dΩ_GW / d ln k  (τ = 300)")
    ax[0, 0].set_title("Espectro de GWs al final de la corrida")
    ax[0, 0].legend(fontsize=8, loc="lower center")
    ax[0, 1].set_ylabel("Δ_χ  (τ = 300)")
    ax[0, 1].set_title("Espectro del campo hijo χ al final")
    for j in (0, 1):
        ax[1, j].axhline(1, color=C_TINTA2, lw=0.8)
        ax[1, j].axhspan(0.8, 1.25, color="#eeeeea", zorder=0)
        ax[1, j].set_yscale("log")
        ax[1, j].set_ylim(0.2, 30)
        ax[1, j].set_xlabel("k (unidades de ω★)")
        ax[1, j].set_ylabel("cociente / N = 128")
    fig.tight_layout()
    fig.savefig(FIG / "espectros_finales.png", dpi=130, bbox_inches="tight")

    # 2) Historia: fracción de energía en χ y energía en GWs
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    for et, r in R.items():
        ax[0].semilogy(r["t"], np.maximum(r["f_chi"], 1e-16), label=et, **estilo[et])
        ax[1].semilogy(r["t_gw"], np.where(r["rhoGW"] > 0, r["rhoGW"], np.nan), label=et, **estilo[et])
    ax[0].set_ylim(1e-14, 1)
    ax[0].set_xlabel("τ")
    ax[0].set_ylabel("fracción de la energía en χ")
    ax[0].set_title("Resonancia: energía transferida a χ")
    ax[1].set_xlabel("τ")
    ax[1].set_ylabel("ρ_GW / ρ")
    ax[1].set_title("Energía total en GWs")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "historias.png", dpi=130, bbox_inches="tight")

    # 3) Hoy: h²Ω_GW(f)
    fig, ax = plt.subplots(figsize=(6.2, 4))
    for et, r in R.items():
        ok = r["h2O"] > 0
        ax.loglog(r["f0"][ok], r["h2O"][ok], label=et, **estilo[et])
    ax.set_ylim(ref["h2O_pico"] * 1e-5, max(r["h2O_pico"] for r in R.values()) * 3)
    ax.set_xlabel("f hoy (Hz)")
    ax.set_ylabel(f"h² Ω_GW hoy (g★ = {G_STAR})")
    ax.set_title("Espectro de GWs hoy")
    ax.legend(fontsize=8, loc="lower center")
    fig.tight_layout()
    fig.savefig(FIG / "gws_hoy.png", dpi=130, bbox_inches="tight")

    # Tabla y resumen
    campos = ["N", "kIR", "tau_res", "rhoGW_fin", "k_pico", "gw_pico", "f_pico", "h2O_pico",
              "friedmann_max", "segundos", "memoria_MB"]
    resumen = {et: {c: r[c] for c in campos} for et, r in R.items()}
    for et, r in R.items():
        d = resumen[et]
        d["k_max"] = float(r["k"][-1])
        d["rhoGW_rel"] = r["rhoGW_fin"] / ref["rhoGW_fin"]
        d["gw_pico_rel"] = r["gw_pico"] / ref["gw_pico"]
        d["costo_rel"] = r["segundos"] / ref["segundos"]
        # diferencia típica del espectro de GWs alrededor del pico (factor 3 a cada lado de k_pico ref)
        kk, cc = en_comun(r["k"], r["gw"], ref["k"], ref["gw"])
        banda = (kk > ref["k_pico"] / 3) & (kk < ref["k_pico"] * 3)
        d["gw_banda_pico_desvio_max"] = float(np.max(np.abs(cc[banda] - 1))) if banda.any() else None
        d["gw_banda_pico_desvio_rms"] = float(np.sqrt(np.mean((cc[banda] - 1) ** 2))) if banda.any() else None
        kk, cc = en_comun(r["k_chi"], r["chi"], ref["k_chi"], ref["chi"])
        banda = (kk > 0.5) & (kk < 4)
        d["chi_k0.5a4_desvio_rms"] = float(np.sqrt(np.mean((cc[banda] - 1) ** 2))) if banda.any() else None
    json.dump(resumen, open(DATA / "convergencia_N" / "resumen.json", "w"), indent=2, ensure_ascii=False)

    print(f"{'corrida':22s} {'k_max':>6s} {'τ_res':>6s} {'ρGW/ρ rel':>9s} {'k_pico':>6s} {'pico rel':>8s} "
          f"{'banda rms':>9s} {'banda max':>9s} {'χ rms':>6s} {'f_pico Hz':>9s} {'Friedm.':>8s} {'tiempo':>8s} {'MB':>5s}")
    for et, d in resumen.items():
        fmt = lambda v: f"{v:9.3f}" if v is not None else f"{'—':>9s}"
        print(f"{et:22s} {d['k_max']:6.2f} {d['tau_res']:6.1f} {d['rhoGW_rel']:9.3f} {d['k_pico']:6.2f} "
              f"{d['gw_pico_rel']:8.3f} {fmt(d['gw_banda_pico_desvio_rms'])} {fmt(d['gw_banda_pico_desvio_max'])} "
              f"{d['chi_k0.5a4_desvio_rms'] or 0:6.3f} {d['f_pico']:9.3g} {d['friedmann_max']:8.1e} "
              f"{d['segundos']/60:7.1f}m {d['memoria_MB']:5.0f}")


if __name__ == "__main__":
    main()
