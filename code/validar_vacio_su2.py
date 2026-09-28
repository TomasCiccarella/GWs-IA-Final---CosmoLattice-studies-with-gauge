"""Validación del parche su2_vacio_transversal.patch (lphi4SU2U1, N = 64, kIR = 0.5 en unidades del control).

Compara dos corridas cortas con la misma semilla: la condición inicial por defecto de CosmoLattice
(data/prueba_su2_vacioT/defecto) y la de vacío transversal en U(1) y SU(2) (data/prueba_su2_vacioT/vacioT).

Chequea:
  1. El doblete (inflatón) arranca con exactamente el mismo ruido en las dos.
  2. La ley de Gauss (U(1) y SU(2)) en t = 0 y durante la prueba.
  3. La energía eléctrica y magnética que agrega el parche, contra la cuenta analítica del vacío de dos
     polarizaciones por autoestado de masa (fotón, Z, W, W; masas y componentes leídas del log).
  4. Que la parte agregada sea casi estacionaria (E·a⁴ y B·a⁴): con masas equivocadas el vacío "chapotea"
     entre E y B (así se descubrió la mezcla A–B^n̂; ver paginas/bitacora.html#vacio-su2).

Uso: python3 code/validar_vacio_su2.py  (escribe data/prueba_su2_vacioT/validacion.json)
"""
import json
import re
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parent.parent
DIR = RAIZ / "data" / "prueba_su2_vacioT"

# Parámetros de la prueba (los del .in, unidades de programa del modelo gauge)
N, KIR, KCUT = 64, 0.707107, 5.65685
LAMBDA = 9e-14  # (omega*/f*)^2 = lambda, porque omega* = sqrt(lambda) f*


def leer(v, nombre):
    return np.loadtxt(DIR / v / nombre, comments="#")


def espectro_t0(v, nombre):
    """Primer bloque (t = 0) de un archivo de espectros: la línea de títulos y los bins hasta la línea vacía."""
    filas = []
    for linea in (DIR / v / nombre).read_text().splitlines()[1:]:
        if not linea.strip():
            break
        filas.append([float(x) for x in linea.split()])
    return np.array(filas)


def autoestados(v):
    """Masas y componentes (A, B1, B2, B3) de los autoestados gauge, leídos del log de la corrida."""
    modos = []
    for linea in (DIR / v / "salida.log").read_text().splitlines():
        m = re.search(r"gauge eigenmode \d+: mass\^2 = ([-\d.e+]+), components \(A, B1, B2, B3\):(.*)", linea)
        if m:
            modos.append((float(m.group(1)), [float(x) for x in m.group(2).split()]))
    return modos


def energia_vacio_por_polarizacion(M2):
    """rho_E y rho_B (unidades de programa) de una polarización en su vacío, sumando los modos de la red."""
    n = np.fft.fftfreq(N, 1.0 / N)
    nx, ny, nz = np.meshgrid(n, n, n, indexing="ij")
    k = KIR * np.sqrt(nx**2 + ny**2 + nz**2)
    dx = 2 * np.pi / (N * KIR)
    klat2 = sum(4 * np.sin(KIR * ni * dx / 2) ** 2 / dx**2 for ni in (nx, ny, nz))
    sel = (k > 0) & (k < KCUT)
    omega = np.sqrt(k[sel] ** 2 + M2)
    V = (N * dx) ** 3
    return LAMBDA * np.sum(omega / 4) / V, LAMBDA * np.sum(klat2[sel] / (4 * omega)) / V


def main():
    res = {}

    # 1. Doblete idéntico en t = 0
    d0 = [leer(v, "average_SU2Doublet_0_0.txt")[0] for v in ("defecto", "vacioT")]
    s0 = [espectro_t0(v, "spectra_norm_SU2Doublet_scalar_0.txt") for v in ("defecto", "vacioT")]
    res["doblete_max_dif_relativa_t0"] = float(
        max(np.max(np.abs(d0[1] - d0[0]) / np.maximum(np.abs(d0[0]), 1e-300)),
            np.max(np.abs(s0[1][:, 1:3] - s0[0][:, 1:3]) / np.maximum(np.abs(s0[0][:, 1:3]), 1e-300))))

    # 2. Gauss
    for g in ("SU2", "U1"):
        for v in ("defecto", "vacioT"):
            a = leer(v, f"average_gauss_{g}_0.txt")
            res[f"gauss_{g}_{v}_t0"] = float(a[0, 1])
            res[f"gauss_{g}_{v}_max"] = float(a[:, 1].max())

    # 3. Energía agregada en t = 0 (columnas de average_energies.txt)
    col = {"EkU1": 7, "EgU1": 8, "EkSU2": 9, "EgSU2": 10}
    e = {v: leer(v, "average_energies.txt") for v in ("defecto", "vacioT")}
    nt = min(len(e["defecto"]), len(e["vacioT"]))
    t = e["defecto"][:nt, 0]
    agregado = {c: e["vacioT"][:nt, i] - e["defecto"][:nt, i] for c, i in col.items()}
    # Esperado: cada autoestado k (masa m_k, componente O_Ak en A) aporta dos polarizaciones de vacío,
    # repartidas entre A (peso O_Ak²) y los colores de SU(2) (peso 1 - O_Ak²).
    modos = autoestados("vacioT")
    res["autoestados"] = [{"m2": m2, "componentes_A_B1_B2_B3": c} for m2, c in modos]
    esperado = dict.fromkeys(col, 0.0)
    for m2, c in modos:
        rE, rB = energia_vacio_por_polarizacion(m2)
        pA = c[0] ** 2
        esperado["EkU1"] += 2 * rE * pA
        esperado["EgU1"] += 2 * rB * pA
        esperado["EkSU2"] += 2 * rE * (1 - pA)
        esperado["EgSU2"] += 2 * rB * (1 - pA)
    res["energia_t0"] = {c: {"agregada": float(agregado[c][0]), "esperada": float(esperado[c]),
                             "cociente": float(agregado[c][0] / esperado[c])} for c in col}

    # 4. Evolución
    # Estacionariedad: E·a⁴ y B·a⁴ de la parte agregada, relativos a t = 0. Un vacío con las masas correctas
    # cambia poco (referencia: lphi4U1, 1,0–1,18 en E y 0,87–1,0 en B); con masas equivocadas "chapotea".
    a4 = leer("vacioT", "average_scale_factor.txt")[:nt, 1] ** 4
    res["estacionariedad"] = {c: (agregado[c] * a4 / agregado[c][0]).round(3).tolist() for c in col}

    (DIR / "validacion.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
