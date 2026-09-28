"""Avance de una corrida de CosmoLattice: porcentaje, velocidad y tiempo restante estimado.

Uso: python3 code/progreso_corrida.py data/<corrida>
Lee el .in (dt, tMax), el último "Step N done" de salida.log y la hora de arranque de progreso.log
(o, si no está, la fecha de salida.log). Solo lee; no toca la corrida.
"""
import re
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

d = Path(sys.argv[1])
par = {}
for linea in next(d.glob("*.in")).read_text().splitlines():
    m = re.match(r"\s*(\w+)\s*=\s*(\S+)", linea)
    if m:
        par[m.group(1)] = m.group(2)
total = round(float(par["tMax"]) / float(par["dt"]))

prog = d / "progreso.log"
if not prog.exists() or "arranque" not in prog.read_text():
    print(f"{d.name}: todavía no arrancó (espera a que terminen las semillas).")
    sys.exit()
tok = prog.read_text().split("arranque: ")[1].split("\n")[0].split()  # salida de `date`: ... PM -03 2026
arranque = datetime.strptime(" ".join(tok[:5] + tok[-1:]), "%a %b %d %I:%M:%S %p %Y")
if "terminada" in prog.read_text():
    print(f"{d.name}: {prog.read_text().strip().splitlines()[-1]}")
    sys.exit()

pasos = [int(p) for p in re.findall(r"Step (\d+) done", (d / "salida.log").read_text())]
paso = pasos[-1] if pasos else 0
transcurrido = datetime.now() - arranque
frac = paso / total
if paso > 0:
    restante = transcurrido * (total - paso) / paso
    fin = datetime.now() + restante
    h, r = divmod(int(restante.total_seconds()), 3600)
    print(f"{d.name}: paso {paso} de {total} ({100 * frac:.1f} %), {transcurrido.total_seconds() / paso:.2f} s/paso; "
          f"faltan ~{h} h {r // 60:02d} min (fin estimado {fin:%d/%m %H:%M}).")
else:
    print(f"{d.name}: arrancó a las {arranque:%H:%M}, todavía sin pasos escritos.")
