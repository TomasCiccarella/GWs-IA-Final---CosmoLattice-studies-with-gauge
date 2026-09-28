"""Genera final-project.pdf a partir de final-project.html (mismo contenido, con las figuras).

Necesita WeasyPrint (`pip install weasyprint`) y las bibliotecas de sistema pango y harfbuzz (en Ubuntu vienen
instaladas). Se usó WeasyPrint 70.0 con el Python del sistema (el de miniconda no encontraba pango).

Uso: python3 code/generar_pdf.py
"""
from pathlib import Path

from weasyprint import CSS, HTML

FINAL = Path(__file__).resolve().parent.parent
IMPRESION = CSS(string="""
@page { size: A4; margin: 16mm 14mm 18mm 14mm;
        @bottom-center { content: "Campos de gauge y ondas gravitacionales — " counter(page) " / " counter(pages);
                         font-size: 8pt; color: #888; } }
body { max-width: none; margin: 0; padding: 0; font-size: 10pt; }
h2 { break-after: avoid; } h3 { break-after: avoid; }
figure, tr { break-inside: avoid; }
.tabla table { break-inside: auto; }
figure img { max-width: 100%; }
th, td { font-size: 8.5pt; }
a { color: #a5651e; text-decoration: none; }
""")

if __name__ == "__main__":
    HTML(FINAL / "final-project.html", base_url=str(FINAL)).write_pdf(FINAL / "final-project.pdf", stylesheets=[IMPRESION])
    print("final-project.pdf escrito")
