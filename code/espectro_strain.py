"""Strain característico hoy (h_c contra f) de los tres modelos, junto a las sensibilidades proyectadas de detectores
de ultra-alta frecuencia. Es la misma figura que GW_Strain_Spectrum.pdf de la tesis de T. Ciccarella (notebook
cosmolattice/GW_Exp.ipynb), con los espectros de este proyecto en lugar de los de la tesis.

Espectros: los de code/produccion_gws.py (semilla 1234, N = 64, kIR = 0.5), llevados a hoy con la misma cuenta
(produccion_gws.gws_hoy: radiación y conservación de la entropía, g★ = 100), y pasados a strain con
    Ω_GW(f) = 2π² f² h_c² / (3 H0²)   →   h_c = √(3/2) · (H0/h) · √(h²Ω_GW) / (π f),
con H0/h = 100 km/s/Mpc = 3.2408×10⁻¹⁸ s⁻¹ (así h_c no depende de h). La parte k > 4 (no convergida) va punteada.

Detectores: los rectángulos y la línea de haces gaussianos son los del notebook de la tesis (mismas coordenadas en
log10 f y log10 h_c y mismos colores), adaptados de Mohanty, Panda y Vidyarthi, arXiv:2503.06858, fig. de
sensibilidades. La curva de cavidades resonantes (EMRC) es la traza aproximada que dejó el mismo notebook (una recta
de (5, −24.75) a (7.75, −33) y una parábola hasta log f = 9); en la figura de la tesis la curva real, con picos,
se dibujó aparte. No se usa para ninguna conclusión: solo ubica los espectros frente a los detectores.

Uso: python3 code/espectro_strain.py   (figura en figures/espectro_strain/, números en data/espectro_strain.json)
"""
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba

FINAL = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(FINAL / "code"))
from produccion_gws import MODELOS, cargar, gws_hoy  # noqa: E402

FIG = FINAL / "figures" / "espectro_strain"
H0_SOBRE_H = 3.2408e-18     # s⁻¹, 100 km/s/Mpc
K_CONV = 4                  # parte confiable del espectro (bitácora 3.3)

# (etiqueta, (log f mín, log h_c mín), ancho, alto, relleno, borde, alpha del relleno, zorder) — de GW_Exp.ipynb
DETECTORES = [
    ("Leviated Sensors", (4.25, -19.25), 1.0, 5, "#4e0011", "#3a000d", 0.4, 1),
    ("ADMX", (4.25, -21), 1.0, 7, "#f1a4be", "#704E5A", 0.2, 0),
    ("BAWs", (6, -20), 2.75, 6, "#4e5363", "#2c2e36", 0.2, 0),
    ("QuEST", (7, -15.75), 0.75, 1, "#0e36b9", "#071c5f", 0.2, 1),
    ("SQMS", (8.75, -23), 0.25, 10, "#0f8d04", "#063d01", 0.2, 0),
    ("ARCADE 2", (9.25, -24), 1.5, 10, "#048d8d", "#013b3d", 0.2, 0),
    ("IAXO SPD", (10.75, -25), 0.25, 15, "#e42020", "#8d0a0a", 0.2, 0),
    ("JURA", (14.25, -31.5), 0.25, 20, "#92e770", "#396934", 0.2, 0),
    ("ALPS II", (14.25, -29), 0.25, 20, "#385791", "#0F2A44", 0.2, 1),
    ("ALPS", (14.25, -25), 0.75, 20, "#9B0CC7", "#310244", 0.2, 2),
]


def strain(f, h2omega):
    return np.sqrt(1.5) * H0_SOBRE_H * np.sqrt(h2omega) / (np.pi * f)


def main():
    mpl.rc("font", family="serif", size=14)
    mpl.rcParams["mathtext.fontset"] = "stix"
    FIG.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.9)
    for nombre, xy, ancho, alto, relleno, borde, alfa, z in DETECTORES:
        ax.add_patch(patches.Rectangle(xy, ancho, alto, facecolor=to_rgba(relleno, alfa),
                                       edgecolor=to_rgba(borde, 1.0), linewidth=2, label=nombre, zorder=z))
    ax.vlines(10.875, -31, -15, color="black", linewidth=1.5, label="EM Gaussian Beams", zorder=2)
    x1 = np.linspace(5, 7.75, 200)
    x2 = np.linspace(7.75, 9, 200)
    ax.plot(np.r_[x1, x2], np.r_[(-33 + 24.75) / (7.75 - 5) * (x1 - 5) - 24.75, 14 / 5 * x2 ** 2 - 461 / 10 * x2 + 1561 / 10],
            color="#314B2E", lw=1.5, label="EMRC (traza aprox.)", zorder=2)

    res = {}
    for m, p in MODELOS.items():
        r = cargar(p)
        f0, h2 = gws_hoy(r)
        k = r["k"]
        ok = (k > 0) & (h2 > 0)
        hc = np.full_like(f0, np.nan)
        hc[ok] = strain(f0[ok], h2[ok])
        conf = ok & (k <= K_CONV + 1e-3)
        resto = ok & (k >= K_CONV - 1e-3)
        ax.plot(np.log10(f0[conf]), np.log10(hc[conf]), color=p["color"], lw=2.2, label=m, zorder=5)
        ax.plot(np.log10(f0[resto]), np.log10(hc[resto]), color=p["color"], lw=1.3, ls=":", zorder=5)
        j = int(np.nanargmax(np.where(conf, hc, np.nan)))
        res[m] = dict(f_hc_max_k_menor_4=float(f0[j]), hc_max_k_menor_4=float(hc[j]),
                      log10_f_rango_confiable=[float(np.log10(f0[conf]).min()), float(np.log10(f0[conf]).max())])

    ax.plot([], [], color="0.4", lw=1.3, ls=":", label="f > 1,3×10⁸ Hz (k > 4, no convergido)")
    ax.set_xlim(4, 16)
    ax.set_ylim(-40, -15)
    ax.set_xticks(np.arange(4, 18, 2))
    ax.set_yticks(np.arange(-40, -14, 3))
    ax.set_xlabel(r"$\log ~ f$  [Hz]", fontsize=18)
    ax.set_ylabel(r"$\log ~ h_c$", fontsize=18)
    ax.legend(loc=(1.02, 0.0), fontsize=11.5)
    fig.tight_layout()
    fig.savefig(FIG / "strain_vs_detectores.png", dpi=150, bbox_inches="tight")

    salida = FINAL / "data" / "espectro_strain.json"
    salida.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
