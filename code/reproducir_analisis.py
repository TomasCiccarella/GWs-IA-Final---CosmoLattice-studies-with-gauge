"""Rehace todos los análisis y verificaciones del proyecto a partir de los datos del repositorio (sin simular).

Corre, en orden, cada script de code/ que produce números, figuras o verificaciones, y al final el control de
procedencia (.claude/hooks/provenance_gate.py). Se detiene en el primer error. Después de correrlo,
`git status` debería mostrar el árbol sin cambios: las figuras y los .json regenerados son los mismos que los
commiteados (salvo diferencias de bytes en los PNG si cambia la versión de matplotlib).

No incluye el notebook del piloto (code/analisis_corridas.py, que es marimo): se abre con
`marimo edit code/analisis_corridas.py` o se corre con `CORRIDA=<carpeta> python3 code/analisis_corridas.py`.

Uso: python3 code/reproducir_analisis.py
"""
import subprocess
import sys
import time
from pathlib import Path

FINAL = Path(__file__).resolve().parent.parent

PASOS = [
    ("Verificaciones simbólicas de unidades, Gauss y cargas (sympy)", ["code/verificar_variables_de_programa.py"]),
    ("Parámetros comparables y Floquet", ["code/analisis_parametros.py", "--tabla"]),
    ("Convergencia en N del control", ["code/convergencia_N.py"]),
    ("Control contra U(1)", ["code/comparar_modelos.py"]),
    ("Crecimiento por modo: U(1) con la condición de CosmoLattice",
     ["code/crecimiento_por_modo_U1.py", "lphi4U1_N64_kIR0.5_VV2"]),
    ("Crecimiento por modo: U(1) con vacío transversal", ["code/crecimiento_por_modo_U1.py", "lphi4U1_vacioT_N64_kIR0.5_VV2"]),
    ("Crecimiento por modo: Z + fotón en SU(2)×U(1)",
     ["code/crecimiento_por_modo_U1.py", "lphi4SU2U1_vacioT_N64_kIR0.5_VV2"]),
    ("Validación del parche SU(2)", ["code/validar_vacio_su2.py"]),
    ("Tres semillas de control y U(1)", ["code/semillas_U1.py"]),
    ("SU(2)×U(1) contra control y U(1)", ["code/analisis_SU2U1.py"]),
    ("Cómo y cuándo se producen las GWs; espectro de hoy", ["code/produccion_gws.py"]),
    ("Forma del espectro: pendientes y pico", ["code/pendientes_espectro.py"]),
]


def main():
    for i, (que, cmd) in enumerate(PASOS, 1):
        t0 = time.time()
        print(f"[{i}/{len(PASOS)}] {que}: {' '.join(cmd)}", flush=True)
        r = subprocess.run([sys.executable, *cmd], cwd=FINAL, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-2000:], r.stderr[-2000:], sep="\n")
            sys.exit(f"Falló el paso {i} ({que}).")
        print(f"      ok ({time.time() - t0:.0f} s)")
    r = subprocess.run([sys.executable, ".claude/hooks/provenance_gate.py"], cwd=FINAL, input="{}",
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr)
        sys.exit("El control de procedencia encontró resultados sin registrar.")
    print("Control de procedencia: todo registrado.")


if __name__ == "__main__":
    main()
