# Para un agente que clona este repositorio

1. Leé `README.md` entero: la pregunta, el enfoque y **"Cómo reproducirlo"** (clonar CosmoLattice en el
   commit fijado, aplicar los parches de `code/parches/` y compilar los builds).
2. Entorno de Python: `pip install -r requirements.txt` (Python 3.12).
3. Todo resultado tiene su registro en `provenance/`: `claims.yaml` (qué se afirma, cómo se verificó y qué lo
   respalda) y `numbers.json` (cada número, con el script que lo produce). El formato está en
   `.claude/provenance/*.md`. Si producís una figura o un número nuevo, registralo ahí: el hook
   `.claude/hooks/provenance_gate.py` no deja terminar un turno sin eso.
4. Los análisis se rehacen desde los datos del repo, sin correr simulaciones: `python3 code/reproducir_analisis.py`
   corre todos los scripts de análisis y de verificación en orden y termina con el control de procedencia.
5. Las simulaciones son caras (control ~11 min, U(1) ~1 h 10 min, SU(2)×U(1) ~9 h 20 min con 8 núcleos y N = 64);
   cada carpeta de `data/` tiene su `.in` y un `LEEME.md` o un `correr.sh` con cómo se lanzó.
6. Idioma: todo lo del proyecto está en español y pensado para que lo entienda alguien que no es de física
   (ver `paginas/bitacora.html`). Los mensajes de commit también van en español.
