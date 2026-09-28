"""¿Cuánto del resultado U(1) vs control es azar? Tres semillas del control y de U(1) con vacío transversal.

Misma red que code/comparar_modelos.py (N = 64, kIR = 0.5 en unidades del control, VV2):
  semilla 1234:       data/convergencia_N/lphi4_N64_kIR0.5_VV2 y data/lphi4U1_vacioT_N64_kIR0.5_VV2
  semillas 2345, 3456: data/semillas/<corrida>_s<semilla>  (lanzadas con data/semillas/correr.sh)

Para cada semilla: τ al 1 % y al 10 % de la energía en el campo hijo/gauge, fracción máxima, ρ_GW/ρ final
y el cociente bin a bin del espectro final de GWs (U(1) / control) en la parte convergida, k = 0.5-2.5.
El cociente se calcula dentro de la misma semilla y también entre promedios de las tres.

Uso: python3 code/semillas_U1.py  (figura en figures/semillas/, números en data/semillas/resumen.json)
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from comparar_modelos import C_AZUL, S2, cargar  # noqa: E402
from convergencia_N import C_AQUA, C_TINTA2, primer_cruce, terminada  # noqa: E402

DATA = FINAL / "data"
FIG = FINAL / "figures" / "semillas"
K_MIN, K_MAX = 0.5, 2.5
ESTILOS = {1234: "-", 2345: "--", 3456: ":"}

MODELOS = {
    "control": dict(carpetas={1234: "convergencia_N/lphi4_N64_kIR0.5_VV2",
                              2345: "semillas/lphi4_N64_kIR0.5_VV2_s2345",
                              3456: "semillas/lphi4_N64_kIR0.5_VV2_s3456"},
                    s=1.0, cols=["E^kin_scal1", "E^grad_scal1"], color=C_AZUL, etiqueta="control λφ⁴ + χ"),
    "U1_vT": dict(carpetas={1234: "lphi4U1_vacioT_N64_kIR0.5_VV2",
                            2345: "semillas/lphi4U1_vacioT_N64_kIR0.5_VV2_s2345",
                            3456: "semillas/lphi4U1_vacioT_N64_kIR0.5_VV2_s3456"},
                  s=S2, cols=["E^kin_U10", "E^grad_U10"], color=C_AQUA, etiqueta="U(1), vacío transversal"),
}


def resumen(r):
    return dict(tau_1pc=primer_cruce(r["t"], r["frac"], 0.01), tau_10pc=primer_cruce(r["t"], r["frac"], 0.1),
                frac_max=float(r["frac"].max()), rhoGW_fin=float(r["rhoGW"][-1]))


def media_desv(x):
    x = np.asarray(x, dtype=float)
    return dict(media=float(x.mean()), desv=float(x.std(ddof=1)) if len(x) > 1 else float("nan"),
                min=float(x.min()), max=float(x.max()))


def main():
    R = {m: {} for m in MODELOS}
    for m, p in MODELOS.items():
        for sem, c in p["carpetas"].items():
            if terminada(DATA / c):
                R[m][sem] = cargar(c, p["s"], p["cols"])
    semillas = sorted(set(R["control"]) & set(R["U1_vT"]))
    print("semillas con las dos corridas terminadas:", semillas)

    k = R["control"][semillas[0]]["k"]
    sel = (k >= K_MIN) & (k <= K_MAX)
    res = {"semillas": semillas, "por_semilla": {}, "entre_semillas": {}}
    cocientes = {}
    for sem in semillas:
        c, u = R["control"][sem], R["U1_vT"][sem]
        assert np.allclose(c["k"], u["k"]), "las grillas de k no coinciden"
        q = u["gw"][sel] / c["gw"][sel]
        cocientes[sem] = q
        res["por_semilla"][sem] = dict(control=resumen(c), U1_vT=resumen(u),
                                       cociente_gw_k05a25=[float(q.min()), float(q.max())],
                                       cociente_gw_k05a25_medio=float(np.exp(np.log(q).mean())),
                                       cociente_rhoGW=float(u["rhoGW"][-1] / c["rhoGW"][-1]))

    # Variación entre semillas de cada modelo, y cociente entre los espectros promedio
    for m in MODELOS:
        res["entre_semillas"][m] = {x: media_desv([resumen(R[m][s])[x] for s in semillas])
                                    for x in ("tau_1pc", "tau_10pc", "frac_max", "rhoGW_fin")}
        g = np.array([R[m][s]["gw"][sel] for s in semillas])
        res["entre_semillas"][m]["dispersion_gw_k05a25"] = float(np.max(g.std(axis=0, ddof=1) / g.mean(axis=0)))
    gc = np.mean([R["control"][s]["gw"][sel] for s in semillas], axis=0)
    gu = np.mean([R["U1_vT"][s]["gw"][sel] for s in semillas], axis=0)
    res["cociente_de_promedios_k05a25"] = [float((gu / gc).min()), float((gu / gc).max())]
    res["cociente_de_promedios_k05a25_medio"] = float(np.exp(np.log(gu / gc).mean()))
    res["cociente_medio_por_semilla"] = media_desv([res["por_semilla"][s]["cociente_gw_k05a25_medio"]
                                                     for s in semillas])
    json.dump(res, open(DATA / "semillas" / "resumen.json", "w"), indent=2, ensure_ascii=False)

    # Figura: espectros finales (izq.) y cociente U(1)/control por semilla (der.)
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for m, p in MODELOS.items():
        for s in semillas:
            r = R[m][s]
            ax[0].loglog(r["k"], np.where(r["gw"] > 0, r["gw"], np.nan), color=p["color"], ls=ESTILOS[s],
                         lw=1.5, label=f"{p['etiqueta']}, semilla {s}")
    for s in semillas:
        ax[1].plot(k[sel], cocientes[s], color=C_TINTA2, ls=ESTILOS[s], lw=1.5, label=f"semilla {s}")
    ax[1].plot(k[sel], gu / gc, color=C_AQUA, lw=2.5, label="cociente de los promedios")
    ax[1].axhline(1, color=C_TINTA2, lw=0.8)
    ax[0].axvspan(4, 30, color="#eeeeea", zorder=0)
    ax[0].text(0.97, 0.05, "k > 4: no convergido", transform=ax[0].transAxes, ha="right", color=C_TINTA2,
               fontsize=8)
    ax[0].set_xlabel("k (unidades del control)")
    ax[0].set_ylabel("dΩ_GW / d ln k (final)")
    ax[0].set_title("Espectro de GWs al final (τ = 300)")
    ax[1].set_xlabel("k (unidades del control)")
    ax[1].set_ylabel("U(1) vacío transversal / control")
    ax[1].set_title(f"Cociente de espectros, k = {K_MIN}–{K_MAX}")
    ax[1].set_ylim(0, None)
    ax[0].legend(fontsize=7)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "U1_vs_control_semillas.png", dpi=130, bbox_inches="tight")

    for s in semillas:
        d = res["por_semilla"][s]
        print(f"semilla {s}: τ1% control {d['control']['tau_1pc']:.1f} U1 {d['U1_vT']['tau_1pc']:.1f} | "
              f"cociente GW k 0.5–2.5 {d['cociente_gw_k05a25'][0]:.2f}–{d['cociente_gw_k05a25'][1]:.2f} "
              f"(medio {d['cociente_gw_k05a25_medio']:.2f}) | ρGW U1/control {d['cociente_rhoGW']:.2f}")
    print("cociente de promedios:", np.round(res["cociente_de_promedios_k05a25"], 2),
          "medio", round(res["cociente_de_promedios_k05a25_medio"], 2))
    for m in MODELOS:
        e = res["entre_semillas"][m]
        print(f"{m}: τ1% {e['tau_1pc']['media']:.1f} ± {e['tau_1pc']['desv']:.1f}, "
              f"ρGW {e['rhoGW_fin']['media']:.2e} ± {e['rhoGW_fin']['desv']:.1e}, "
              f"dispersión máx. del espectro k 0.5–2.5 entre semillas {100 * e['dispersion_gw_k05a25']:.0f} %")


if __name__ == "__main__":
    main()
