"""Genera los dos PDF del proyecto a partir de sus páginas HTML (mismo contenido, con las figuras).

  informe-definitivo.pdf  <- informe-definitivo.html  (el informe breve pedido: 5 páginas como máximo)
  informe-extendido.pdf   <- final-project.html       (la página completa del proyecto)

Necesita WeasyPrint (`pip install weasyprint`) y las bibliotecas de sistema pango y harfbuzz (en Ubuntu vienen
instaladas). Se usó WeasyPrint 70.0 con el Python del sistema (el de miniconda no encontraba pango).
Verifica que el informe definitivo no pase de 5 páginas.

Uso: python3 code/generar_pdf.py
"""
import sys
from pathlib import Path

from weasyprint import CSS, HTML

FINAL = Path(__file__).resolve().parent.parent
PIE = "Campos de gauge y ondas gravitacionales"


def estilo(tamano, pie):
    return CSS(string=f"""
@page {{ size: A4; margin: 14mm 13mm 16mm 13mm;
        @bottom-center {{ content: "{pie} — " counter(page) " / " counter(pages); font-size: 8pt; color: #888; }} }}
body {{ max-width: none; margin: 0; padding: 0; font-size: {tamano}; }}
h2, h3 {{ break-after: avoid; }}
figure, tr {{ break-inside: avoid; }}
.tabla table {{ break-inside: auto; }}
figure img {{ max-width: 100%; }}
a {{ color: #a5651e; text-decoration: none; }}
""")


PDFS = [
    ("informe-definitivo.html", "informe-definitivo.pdf", estilo("9.2pt", PIE + " · informe definitivo"), 5),
    ("final-project.html", "informe-extendido.pdf", estilo("10pt", PIE + " · informe extendido"), None),
]

if __name__ == "__main__":
    for html, pdf, css, max_paginas in PDFS:
        doc = HTML(FINAL / html, base_url=str(FINAL)).render(stylesheets=[css])
        doc.write_pdf(FINAL / pdf)
        n = len(doc.pages)
        print(f"{pdf}: {n} páginas")
        if max_paginas and n > max_paginas:
            sys.exit(f"{pdf} tiene {n} páginas; el máximo es {max_paginas}.")
