"""Forma del espectro final de GWs: pendientes a cada lado del pico, posición y altura del pico, en las 7 corridas.

Responde la parte de la pregunta del proyecto sobre la FORMA del espectro (bitácora §3.9). Corridas en la misma
red (N = 64, kIR = 0.5 en unidades del control): control y U(1) con vacío transversal con tres semillas (1234,
2345, 3456; data/semillas/) y SU(2)×U(1) con una (1234).

Ajuste dΩ_GW/d ln k ∝ k^n por mínimos cuadrados en log-log:
  lado IR (debajo del pico):    k = 0.5, 1.0, 1.5  (los tres primeros bins; la caja no llega más abajo)
  lado UV cercano (sobre el pico): k = 2.0 a 4.0  (hasta donde el espectro está convergido, bitácora §3.3)
Con tres puntos el ajuste IR es grueso; la dispersión entre semillas da la barra de error honesta.

Uso: python3 code/pendientes_espectro.py   (números en data/pendientes_espectro.json, figura en
     figures/pendientes_espectro.png)
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from comparar_modelos import C_AZUL, S2  # noqa: E402
from convergencia_N import C_AQUA, C_NARANJA, C_TINTA2, leer_espectros  # noqa: E402

DATA = FINAL / "data"
IR = (0.5, 1.5)
UV = (2.0, 4.0)
CORRIDAS = {
    "control": (1.0, C_AZUL, {1234: "convergencia_N/lphi4_N64_kIR0.5_VV2",
                              2345: "semillas/lphi4_N64_kIR0.5_VV2_s2345",
                              3456: "semillas/lphi4_N64_kIR0.5_VV2_s3456"}),
    "U(1)": (S2, C_AQUA, {1234: "lphi4U1_vacioT_N64_kIR0.5_VV2",
                          2345: "semillas/lphi4U1_vacioT_N64_kIR0.5_VV2_s2345",
                          3456: "semillas/lphi4U1_vacioT_N64_kIR0.5_VV2_s3456"}),
    "SU(2)×U(1)": (S2, C_NARANJA, {1234: "lphi4SU2U1_vacioT_N64_kIR0.5_VV2"}),
}


def pendiente(k, y, rango):
    m = (k >= rango[0] - 1e-3) & (k <= rango[1] + 1e-3) & (y > 0)
    n, _ = np.polyfit(np.log(k[m]), np.log(y[m]), 1)
    return float(n)


def main():
    res = {}
    fig, ax = plt.subplots(figsize=(7, 4.2))
    estilos = {1234: "-", 2345: "--", 3456: ":"}
    for m, (s, color, carpetas) in CORRIDAS.items():
        res[m] = {"por_semilla": {}}
        for sem, c in carpetas.items():
            esp = leer_espectros(DATA / c / "spectra_energy_gws.txt")
            k, y = esp[-1, :, 0] / s, esp[-1, :, 1]
            sel = (k > 0) & (k <= 4 + 1e-3)
            jp = int(np.argmax(np.where(sel, y, 0)))
            res[m]["por_semilla"][sem] = dict(n_IR=pendiente(k, y, IR), n_UV=pendiente(k, y, UV),
                                              k_pico=float(k[jp]), altura_pico=float(y[jp]))
            ax.loglog(k[sel], y[sel], color=color, ls=estilos[sem], lw=1.5, label=f"{m}, semilla {sem}")
        for x in ("n_IR", "n_UV", "k_pico", "altura_pico"):
            v = np.array([d[x] for d in res[m]["por_semilla"].values()])
            res[m][x] = dict(media=float(v.mean()), desv=float(v.std(ddof=1)) if len(v) > 1 else None,
                             min=float(v.min()), max=float(v.max()))
    for (a, b) in ((0.5, 1.5), (2.0, 4.0)):
        ax.axvspan(a, b, color="#eeeeea", zorder=0)
    ax.text(0.87, 0.03, "ajuste IR", color=C_TINTA2, fontsize=8, ha="center", transform=ax.get_xaxis_transform())
    ax.text(2.83, 0.03, "ajuste UV", color=C_TINTA2, fontsize=8, ha="center", transform=ax.get_xaxis_transform())
    ax.set_xlabel("k (unidades del control)")
    ax.set_ylabel("dΩ_GW / d ln k (final)")
    ax.set_title("Espectro final en la parte convergida (k ≤ 4)")
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(FINAL / "figures" / "pendientes_espectro.png", dpi=130, bbox_inches="tight")
    json.dump(res, open(DATA / "pendientes_espectro.json", "w"), indent=2, ensure_ascii=False)
    for m, d in res.items():
        f = lambda x: f"{d[x]['media']:.2f}" + (f" ± {d[x]['desv']:.2f}" if d[x]["desv"] is not None else "")  # noqa: E731
        print(f"{m:11s} n_IR {f('n_IR')}  n_UV {f('n_UV')}  k_pico {f('k_pico')}  "
              f"altura {d['altura_pico']['media']:.2e} ({d['altura_pico']['min']:.2e}-{d['altura_pico']['max']:.2e})")
        print("   ", {s: {x: round(v, 2) if x != 'altura_pico' else f"{v:.2e}" for x, v in e.items()}
                    for s, e in d["por_semilla"].items()})


if __name__ == "__main__":
    main()
