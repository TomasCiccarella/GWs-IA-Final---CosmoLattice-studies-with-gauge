"""SU(2)×U(1) con vacío transversal contra el control y contra U(1), en la misma red y con la misma semilla.

Corridas (N = 64, kIR = 0.5 en unidades del control, VV2, semilla 1234):
  control:          data/convergencia_N/lphi4_N64_kIR0.5_VV2
  U(1):             data/lphi4U1_vacioT_N64_kIR0.5_VV2      (parche u1_vacio_transversal)
  SU(2)×U(1):       data/lphi4SU2U1_vacioT_N64_kIR0.5_VV2   (parches u1_ y su2_vacio_transversal)

Todo en unidades del control: k = k_gauge/√2, τ = √2 t_gauge.

Separación de autoestados. Con q_A = q_B = 60, alrededor del doblete homogéneo los autoestados de masa son
fotón γ = (A + B¹)/√2 (m² = 0), Z = (A − B¹)/√2 (m² = 240 gauge, q = 120 del control) y W = B², B³
(m² = 120, q = 60); ver bitacora.html#vacio-su2. Entonces A = (Z + γ)/√2 y B¹ = (γ − Z)/√2, y
    ρ_U1 = (ρ_Z + ρ_γ)/2 + (cruzado Z·γ),     ρ_SU2 = (ρ_Z + ρ_γ)/2 − (cruzado Z·γ) + ρ_W.
Así ρ_Z + ρ_γ ≈ 2 ρ_U1 y ρ_W ≈ ρ_SU2 − ρ_U1, exactos salvo el término cruzado (que en promedio se anula si
Z y γ no están correlacionados) y válidos mientras la base de autoestados tenga sentido (fondo homogéneo).

Cuidado con la energía magnética de SU(2): se mide con 2 − tr(plaqueta) y tiene un piso de redondeo mientras
los campos son de vacío (ver data/prueba_su2_vacioT/LEEME.md); por eso los ajustes de crecimiento usan
ventanas donde la señal está muy por encima de ese piso, y se comparan eléctrico y magnético.

Uso: python3 code/analisis_SU2U1.py   (figuras en figures/lphi4SU2U1_vacioT_N64_kIR0.5_VV2/,
                                        números en data/lphi4SU2U1_vacioT_N64_kIR0.5_VV2/analisis.json)
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from analisis_parametros import floquet_mu  # noqa: E402
from convergencia_N import C_AQUA, C_NARANJA, C_TINTA2, leer_espectros, leer_promedios, primer_cruce  # noqa: E402

DATA = FINAL / "data"
SU2 = DATA / "lphi4SU2U1_vacioT_N64_kIR0.5_VV2"
U1 = DATA / "lphi4U1_vacioT_N64_kIR0.5_VV2"
CONTROL = DATA / "convergencia_N" / "lphi4_N64_kIR0.5_VV2"
FIG = FINAL / "figures" / SU2.name
S2 = np.sqrt(2)
A_CONF = 1.138          # amplitud conforme del inflatón medida en el piloto (bitácora 3.1)
C_AZUL = "#2a78d6"
K_MIN, K_MAX = 0.5, 2.5  # parte convergida del espectro de GWs (bitácora 3.3)
VENTANAS_LINEALES = [(20, 40), (40, 60), (30, 70)]
VENTANA_W_TARDIA = (100, 125)


def cargar(D, s):
    en = leer_promedios(D / "average_energies.txt")
    gw = leer_promedios(D / "average_energies_gws.txt")
    esp = leer_espectros(D / "spectra_energy_gws.txt")
    t_esp = np.atleast_1d(np.loadtxt(D / "average_spectra_times.txt"))[: len(esp)] * s
    return dict(en=en, t=en["t"] * s, T=en["E_tot"], t_gw=gw["t"] * s, rhoGW=gw["rhoGW_over_rho"],
                t_esp=t_esp, k=esp[-1, :, 0] / s, gw=esp[-1, :, 1], esp=esp[:, :, 1])


def mu_energia(t, rho, ventana):
    """μ de la componente que más crece: ρ ∝ e^{2μτ}, así que μ = ½ d ln ρ / dτ (ajuste lineal en la ventana)."""
    m = (t >= ventana[0]) & (t <= ventana[1]) & (rho > 0)
    return float(np.polyfit(t[m], 0.5 * np.log(rho[m]), 1)[0]) if m.sum() > 3 else float("nan")


def floquet_banda(q):
    K = np.linspace(0.02, 4, 160)
    mu = np.array([A_CONF * floquet_mu(k / A_CONF, q, pasos=1500) for k in K])
    i = int(np.argmax(mu))
    dentro = K[mu > mu[i] / 2]
    return dict(mu_max=float(mu[i]), k_max=float(K[i]), k_media_altura=[float(dentro.min()), float(dentro.max())],
                K=K, mu=mu)


def main():
    C, U, S = cargar(CONTROL, 1.0), cargar(U1, S2), cargar(SU2, S2)
    e = S["en"]
    # Fracciones de la energía total (cinética + gradiente de los campos hijos o gauge)
    f = {
        "control: χ": (C["t"], (C["en"]["E^kin_scal1"] + C["en"]["E^grad_scal1"]) / C["T"]),
        "U(1): A": (U["t"], (U["en"]["E^kin_U10"] + U["en"]["E^grad_U10"]) / U["T"]),
        "SU(2)×U(1): todos los gauge": (S["t"], (e["E^kin_U10"] + e["E^grad_U10"] + e["E^kin_SU2"] + e["E^grad_SU2"]) / S["T"]),
        "SU(2)×U(1): Z + fotón (2 ρ_U1)": (S["t"], 2 * (e["E^kin_U10"] + e["E^grad_U10"]) / S["T"]),
        "SU(2)×U(1): W (ρ_SU2 − ρ_U1)": (S["t"], (e["E^kin_SU2"] + e["E^grad_SU2"] - e["E^kin_U10"] - e["E^grad_U10"]) / S["T"]),
    }
    res = {"fracciones": {}}
    for et, (t, y) in f.items():
        res["fracciones"][et] = dict(tau_1pc=primer_cruce(t, y, 0.01), tau_10pc=primer_cruce(t, y, 0.1),
                                     max=float(y.max()), tau_max=float(t[np.argmax(y)]), final=float(y[-1]))
    res["fracciones"]["SU(2)×U(1): gradiente del doblete (final)"] = float(e["E^grad_SU2matter"][-1] / S["T"][-1])
    res["fracciones"]["SU(2)×U(1): U(1) eléctrico/magnético final"] = [float(e["E^kin_U10"][-1] / S["T"][-1]),
                                                                       float(e["E^grad_U10"][-1] / S["T"][-1])]

    # Crecimiento: μ de la energía magnética (transversal pura). W con el magnético y con el eléctrico.
    BZ = 2 * e["E^grad_U10"]
    BW = e["E^grad_SU2"] - e["E^grad_U10"]
    EW = e["E^kin_SU2"] - e["E^kin_U10"]
    q120, q60 = floquet_banda(120), floquet_banda(60)
    res["floquet"] = {q: {x: v for x, v in b.items() if x not in ("K", "mu")} for q, b in ((120, q120), (60, q60))}
    res["mu"] = {}
    for v in VENTANAS_LINEALES + [VENTANA_W_TARDIA]:
        res["mu"][f"{v[0]}-{v[1]}"] = dict(
            control_chi=mu_energia(C["t"], C["en"]["E^grad_scal1"], v),
            U1_A_magnetico=mu_energia(U["t"], U["en"]["E^grad_U10"], v),
            SU2U1_Zgamma_magnetico=mu_energia(S["t"], BZ, v),
            SU2U1_W_magnetico=mu_energia(S["t"], BW, v),
            SU2U1_W_electrico=mu_energia(S["t"], EW, v))
    # ρ_SU2 − ρ_U1 negativo tiene dos causas distintas: antes de τ = 40, el piso de redondeo del magnético SU(2)
    # (campos de vacío); en la fase no lineal (τ ~ 75-110), el término cruzado Z·γ (Z y fotón correlacionados).
    temprano = S["t"] < 40
    res["piso_B_SU2"] = float(-np.min((BW / S["T"])[temprano]))
    nl = (S["t"] > 60) & (S["t"] < 115)
    res["cruzado_Zgamma"] = dict(min_W_magnetico=float(np.min((BW / S["T"])[nl])),
                                 min_W_electrico=float(np.min((EW / S["T"])[nl])),
                                 tau_min=float(S["t"][nl][np.argmin((BW / S["T"])[nl])]),
                                 rho_U1_en_ese_tau=float(((e["E^kin_U10"] + e["E^grad_U10"]) / S["T"])[nl][np.argmin((BW / S["T"])[nl])]))

    # GWs
    sel = (C["k"] >= K_MIN) & (C["k"] <= K_MAX)
    assert np.allclose(C["k"], S["k"]) and np.allclose(C["k"], U["k"])
    res["gw"] = {}
    for et, r in (("control", C), ("U(1)", U), ("SU(2)×U(1)", S)):
        jp = int(np.argmax(np.where(r["k"] <= 4, r["gw"], 0)))
        fin = r["rhoGW"][-1]
        res["gw"][et] = dict(rhoGW_fin=float(fin), k_pico_k_menor_4=float(r["k"][jp]),
                             gw_pico_k_menor_4=float(r["gw"][jp]),
                             tau_mitad_rhoGW=primer_cruce(r["t_gw"], np.maximum(r["rhoGW"], 1e-300), fin / 2))
    for et, num, den in (("SU2U1/control", S, C), ("SU2U1/U1", S, U), ("U1/control", U, C)):
        q = num["gw"][sel] / den["gw"][sel]
        res["gw"][et] = dict(k=num["k"][sel].tolist(), cociente=np.round(q, 3).tolist(),
                             rango=[float(q.min()), float(q.max())], medio=float(np.exp(np.log(q).mean())),
                             rhoGW=float(num["rhoGW"][-1] / den["rhoGW"][-1]))

    # Calidad numérica
    g2 = np.loadtxt(SU2 / "average_gauss_SU2_0.txt")
    g1 = np.loadtxt(SU2 / "average_gauss_U1_0.txt")
    fr = np.loadtxt(SU2 / "average_energy_conservation.txt")
    res["calidad"] = dict(gauss_SU2_max=float(g2[:, 1].max()), gauss_SU2_final=float(g2[-1, 1]),
                          gauss_U1_max=float(g1[:, 1].max()), friedmann_max=float(np.abs(fr[:, 1]).max()))
    json.dump(res, open(SU2 / "analisis.json", "w"), indent=2, ensure_ascii=False)

    # ---------------- Figuras ----------------
    FIG.mkdir(parents=True, exist_ok=True)
    estilo = {"control: χ": (C_AZUL, "-"), "U(1): A": (C_AQUA, "-"),
              "SU(2)×U(1): todos los gauge": (C_NARANJA, "-"),
              "SU(2)×U(1): Z + fotón (2 ρ_U1)": (C_NARANJA, "--"),
              "SU(2)×U(1): W (ρ_SU2 − ρ_U1)": (C_NARANJA, ":")}

    # 1. Transferencia de energía
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for et, (t, y) in f.items():
        c, ls = estilo[et]
        ax[0].semilogy(t, np.where(y > 0, y, np.nan), color=c, ls=ls, lw=1.6, label=et)
        ax[1].plot(t, y, color=c, ls=ls, lw=1.6, label=et)
    ax[0].axhline(res["piso_B_SU2"], color=C_TINTA2, lw=0.7, ls="-.")
    ax[0].text(150, res["piso_B_SU2"] * 1.6, "piso de redondeo de B_SU2 (τ < 40)", color=C_TINTA2, fontsize=7)
    ax[0].axvspan(75, 110, color="#eeeeea", zorder=0)
    ax[0].text(92, 3e-12, "W no separable\n(término cruzado Z·γ)", ha="center", color=C_TINTA2, fontsize=7)
    ax[0].set_ylim(1e-13, 1.5)
    ax[0].set_title("Fracción de la energía en los campos hijos / gauge")
    ax[1].set_title("Lo mismo, escala lineal")
    for a in ax:
        a.set_xlabel("τ (unidades del control)")
    ax[0].set_ylabel("ρ_campo / ρ_total")
    ax[1].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "transferencia_energia.png", dpi=130, bbox_inches="tight")

    # 2. Crecimiento: ½ ln ρ_B contra τ con las pendientes de Floquet
    fig, ax = plt.subplots(figsize=(7, 4.2))
    series = (("control: χ (gradiente)", C["t"], C["en"]["E^grad_scal1"] / C["T"], C_AZUL, "-"),
              ("U(1): A (magnético)", U["t"], U["en"]["E^grad_U10"] / U["T"], C_AQUA, "-"),
              ("SU(2)×U(1): Z + fotón (magnético)", S["t"], BZ / S["T"], C_NARANJA, "--"),
              ("SU(2)×U(1): W (magnético)", S["t"], BW / S["T"], C_NARANJA, ":"))
    for et, t, y, c, ls in series:
        ax.plot(t, 0.5 * np.log(np.where(y > 0, y, np.nan)), color=c, ls=ls, lw=1.5, label=et)
    for q, b, ls in ((120, q120, "-"), (60, q60, "--")):
        i0 = np.argmin(abs(S["t"] - 20)); t0, y0 = 20, 0.5 * np.log(BZ[i0] / S["T"][i0])
        tt = np.array([t0, t0 + 45])
        ax.plot(tt, y0 + b["mu_max"] * (tt - t0), color=C_TINTA2, ls=ls, lw=0.9,
                label=f"pendiente Floquet q = {q} (μ = {b['mu_max']:.2f})")
    ax.set_xlim(0, 160)
    ax.set_xlabel("τ (unidades del control)")
    ax.set_ylabel("½ ln(ρ_campo / ρ_total)")
    ax.set_title("Crecimiento de la energía (magnética para los gauge)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "crecimiento_energia.png", dpi=130, bbox_inches="tight")

    # 3. Ondas gravitacionales
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for et, r, c in (("control", C, C_AZUL), ("U(1)", U, C_AQUA), ("SU(2)×U(1)", S, C_NARANJA)):
        ax[0].semilogy(r["t_gw"], np.where(r["rhoGW"] > 0, r["rhoGW"], np.nan), color=c, lw=1.6, label=et)
        ax[1].loglog(r["k"], np.where(r["gw"] > 0, r["gw"], np.nan), color=c, lw=1.6, label=et)
    for et, c, ls in (("SU2U1/control", C_NARANJA, "-"), ("SU2U1/U1", C_NARANJA, "--"), ("U1/control", C_AQUA, "-")):
        ax[2].plot(res["gw"][et]["k"], res["gw"][et]["cociente"], color=c, ls=ls, lw=1.6, marker="o", ms=4,
                   label=et.replace("SU2U1", "SU(2)×U(1)").replace("U1", "U(1)"))
    ax[2].axhline(1, color=C_TINTA2, lw=0.8)
    ax[1].axvspan(4, 30, color="#eeeeea", zorder=0)
    ax[1].text(0.97, 0.05, "k > 4: no convergido", transform=ax[1].transAxes, ha="right", color=C_TINTA2, fontsize=8)
    ax[0].set_xlabel("τ (unidades del control)")
    ax[0].set_ylabel("ρ_GW / ρ")
    ax[0].set_title("Energía total en GWs")
    ax[1].set_xlabel("k (unidades del control)")
    ax[1].set_ylabel("dΩ_GW / d ln k (final)")
    ax[1].set_title("Espectro de GWs al final (τ = 300)")
    ax[2].set_xlabel("k (unidades del control)")
    ax[2].set_ylabel("cociente de espectros")
    ax[2].set_title(f"Cocientes en k = {K_MIN}–{K_MAX} (una semilla)")
    ax[2].set_ylim(0, None)
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "ondas_gravitacionales.png", dpi=130, bbox_inches="tight")

    # Resumen en pantalla
    for et, d in res["fracciones"].items():
        if isinstance(d, dict):
            print(f"{et:34s} τ1% {d['tau_1pc']:6.1f}  τ10% {d['tau_10pc']:6.1f}  máx {d['max']:.3f} (τ {d['tau_max']:.0f})  final {d['final']:.3f}")
        else:
            print(f"{et:34s} {d}")
    print("Floquet:", res["floquet"])
    for v, d in res["mu"].items():
        print(f"μ τ {v}: " + ", ".join(f"{x} {y:.3f}" for x, y in d.items()))
    print("piso B_SU2:", res["piso_B_SU2"], " cruzado Z·γ:", res["cruzado_Zgamma"])
    for et, d in res["gw"].items():
        print(et, {x: y for x, y in d.items() if x not in ("k",)})
    print("calidad:", res["calidad"])


if __name__ == "__main__":
    main()
