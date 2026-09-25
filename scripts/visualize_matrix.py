# Copyright 2026 Project Contributors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Interactive PFM / PPM / PWM Matrix Heatmap Explorer and Multi-Format Exporter.

Transforms raw Position Frequency Matrices (PFMs) into normalized Position
Probability Matrices (PPMs) and log-odds Position Weight Matrices (PWMs) with
configurable per-base pseudocounts and background GC frequencies. Generates
structured JSON metrics and interactive diverging-color HTML5 heatmaps with
MEME, TRANSFAC, JASPAR, and YAML export capabilities.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

from jaspar_api import (
    ACADEMIC_CITATION_NOTICE,
    SCIENTIFIC_DISCLAIMER,
    compute_ppm_and_pwm,
    format_motif_output,
    load_mock_fixture,
    verify_workspace_license,
)
from visualize_logo import compute_logo_metrics


def build_matrix_analysis(
    data: Dict[str, Any],
    pseudocount: float = 0.8,
    gc_content: float = 0.50,
) -> Dict[str, Any]:
    """Computes PFM, PPM, PWM, positional entropy/bits, and multi-format exports."""
    pfm = data["pfm"]
    calc = compute_ppm_and_pwm(pfm, pseudocount=pseudocount, gc_content=gc_content)
    logo_info = compute_logo_metrics(pfm, pseudocount=pseudocount, reverse_complement=False)

    entropy_bits = [col["entropy_h_i"] for col in logo_info["columns"]]
    info_content_bits = [col["info_r_i"] for col in logo_info["columns"]]

    exports = {
        "jaspar": format_motif_output(data, "jaspar"),
        "meme": format_motif_output(data, "meme"),
        "transfac": format_motif_output(data, "transfac"),
        "yaml": format_motif_output(data, "yaml"),
        "pfm": format_motif_output(data, "pfm"),
    }

    return {
        "matrix_id": data.get("matrix_id", "MA0139.1"),
        "name": data.get("name", "CTCF"),
        "length": calc["length"],
        "pseudocount": pseudocount,
        "gc_content": gc_content,
        "background": calc["background"],
        "column_totals": calc["column_totals"],
        "pfm": pfm,
        "ppm": calc["ppm"],
        "pwm": calc["pwm"],
        "entropy_bits": entropy_bits,
        "info_content_bits": info_content_bits,
        "total_bits": logo_info["total_bits"],
        "consensus": logo_info["consensus"],
        "exports": exports,
        "citations": ACADEMIC_CITATION_NOTICE,
        "disclaimer": SCIENTIFIC_DISCLAIMER,
    }


def render_matrix_html_widget(analysis: Dict[str, Any]) -> str:
    """Renders an interactive HTML5 Matrix Heatmap Explorer with diverging PWM colors and format exporters."""
    title = f"{analysis['matrix_id']} ({analysis['name']}) Matrix Explorer"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, sans-serif; background: #F8FAFC; color: #0F172A; margin: 24px; }}
    .card {{ background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 20px; max-width: 1180px; margin: 0 auto; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
    .controls {{ display: flex; flex-wrap: wrap; gap: 16px; align-items: center; justify-content: space-between; margin-bottom: 16px; background: #F1F5F9; padding: 12px 16px; border-radius: 8px; }}
    .btn-group {{ display: flex; gap: 8px; flex-wrap: wrap; }}
    .btn {{ background: #0F172A; color: #FFFFFF; border: none; border-radius: 6px; padding: 7px 12px; font-size: 12px; font-weight: 600; cursor: pointer; }}
    .btn.secondary {{ background: #FFFFFF; color: #0F172A; border: 1px solid #CBD5E1; }}
    .table-wrap {{ overflow-x: auto; border: 1px solid #E2E8F0; border-radius: 8px; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 12px; font-family: ui-monospace, monospace; }}
    th, td {{ border: 1px solid #E2E8F0; padding: 6px 8px; text-align: right; }}
    th {{ background: #F8FAFC; color: #334155; font-weight: 700; text-align: center; }}
    .base-label {{ font-weight: 700; text-align: center; background: #F8FAFC; width: 48px; }}
    pre.export-box {{ background: #0F172A; color: #F8FAFC; padding: 14px; border-radius: 8px; font-size: 12px; overflow-x: auto; margin-top: 14px; }}
    .footer {{ font-size: 12px; color: #475569; margin-top: 16px; border-top: 1px solid #E2E8F0; padding-top: 12px; }}
  </style>
</head>
<body>
  <div class="card">
    <h2>{html.escape(title)}</h2>
    <div class="controls">
      <div class="btn-group">
        <button class="btn" onclick="setMatrixMode('pwm')">PWM (Log2-Odds)</button>
        <button class="btn secondary" onclick="setMatrixMode('ppm')">PPM (Probabilities)</button>
        <button class="btn secondary" onclick="setMatrixMode('pfm')">PFM (Raw Counts)</button>
      </div>
      <div>
        <label>Pseudocount: <strong id="pseudoVal">{analysis['pseudocount']}</strong></label>
        <input type="range" id="pseudoSlider" min="0.1" max="5.0" step="0.1" value="{analysis['pseudocount']}" oninput="recalculateMatrix()"/>
        &nbsp;&nbsp;
        <label>Background GC: <strong id="gcVal">{analysis['gc_content']:.2f}</strong></label>
        <input type="range" id="gcSlider" min="0.20" max="0.80" step="0.02" value="{analysis['gc_content']}" oninput="recalculateMatrix()"/>
      </div>
      <div class="btn-group">
        <button class="btn secondary" onclick="showExport('jaspar')">JASPAR</button>
        <button class="btn secondary" onclick="showExport('meme')">MEME</button>
        <button class="btn secondary" onclick="showExport('transfac')">TRANSFAC</button>
        <button class="btn secondary" onclick="showExport('yaml')">YAML</button>
      </div>
    </div>
    <div class="table-wrap" id="matrixTableContainer"></div>
    <pre class="export-box" id="exportBox"></pre>
    <div class="footer">
      <p><strong>Academic Citation:</strong> Ovek Baydar, D. et al. JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models. Nucleic Acids Research 54, D184-D193 (2026). doi:10.1093/nar/gkaf1209. Licensed under CC BY 4.0.</p>
      <p><strong>Scientific Disclaimer:</strong> {html.escape(SCIENTIFIC_DISCLAIMER)}</p>
    </div>
  </div>
  <script>
    const MATRIX_STATE = {json.dumps(analysis)};
    let currentMode = 'pwm';
    let currentExport = 'jaspar';

    function divergingColor(val, mode) {{
      if (mode === 'pfm') return '#FFFFFF';
      if (mode === 'ppm') {{
        const alpha = Math.min(0.85, Math.max(0.05, val));
        return `rgba(59, 130, 246, ${{alpha.toFixed(2)}})`;
      }}
      // PWM diverging color scale: negative -> light red, 0 -> neutral white/gray, positive -> blue
      if (val > 0.05) {{
        const alpha = Math.min(0.85, val / 2.0);
        return `rgba(59, 130, 246, ${{alpha.toFixed(2)}})`;
      }} else if (val < -0.05) {{
        const alpha = Math.min(0.75, Math.abs(val) / 4.0);
        return `rgba(239, 68, 68, ${{alpha.toFixed(2)}})`;
      }}
      return '#F8FAFC';
    }}

    function recalculateMatrix() {{
      const pseudo = parseFloat(document.getElementById('pseudoSlider').value);
      const gc = parseFloat(document.getElementById('gcSlider').value);
      document.getElementById('pseudoVal').textContent = pseudo.toFixed(1);
      document.getElementById('gcVal').textContent = gc.toFixed(2);

      const bg = {{ A: (1 - gc) / 2, T: (1 - gc) / 2, G: gc / 2, C: gc / 2 }};
      const L = MATRIX_STATE.length;
      const bases = ['A', 'C', 'G', 'T'];

      let htmlRows = '<table><thead><tr><th>Base</th>';
      for (let i = 1; i <= L; i++) htmlRows += `<th>${{i}}</th>`;
      htmlRows += '</tr></thead><tbody>';

      for (const b of bases) {{
        htmlRows += `<tr><td class="base-label">${{b}}</td>`;
        for (let i = 0; i < L; i++) {{
          const n_i = MATRIX_STATE.pfm.A[i] + MATRIX_STATE.pfm.C[i] + MATRIX_STATE.pfm.G[i] + MATRIX_STATE.pfm.T[i];
          const raw = MATRIX_STATE.pfm[b][i];
          const prob = (raw + pseudo) / (n_i + 4 * pseudo);
          const weight = Math.log2(prob / bg[b]);
          const displayVal = currentMode === 'pfm' ? raw.toFixed(0) : (currentMode === 'ppm' ? prob.toFixed(3) : weight.toFixed(2));
          const bgColor = divergingColor(currentMode === 'pwm' ? weight : prob, currentMode);
          htmlRows += `<td style="background:${{bgColor}}">${{displayVal}}</td>`;
        }}
        htmlRows += '</tr>';
      }}
      htmlRows += '</tbody></table>';
      document.getElementById('matrixTableContainer').innerHTML = htmlRows;
      document.getElementById('exportBox').textContent = MATRIX_STATE.exports[currentExport];
    }}

    function setMatrixMode(mode) {{ currentMode = mode; recalculateMatrix(); }}
    function showExport(fmt) {{ currentExport = fmt; recalculateMatrix(); }}
    recalculateMatrix();
  </script>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    """Constructs CLI argument parser for visualize_matrix.py."""
    parser = argparse.ArgumentParser(
        prog="visualize_matrix.py",
        description="Transform raw PFM counts into normalized PPMs, PWMs, and interactive heatmaps.",
    )
    parser.add_argument("--matrix-json", help="Path to JSON file with matrix PFM.")
    parser.add_argument("--pseudocount", type=float, default=0.8, help="Per-base pseudocount (default: 0.8).")
    parser.add_argument("--gc-content", type=float, default=0.50, help="Background GC frequency (default: 0.50).")
    parser.add_argument("--output-json", help="Export normalized PPM, PWM, and position entropy to JSON.")
    parser.add_argument("--output-html", help="Export interactive tabular heatmap component.")
    parser.add_argument("--mock", action="store_true", help="Offline mock fallback.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for visualize_matrix.py."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        verify_workspace_license(quiet=True)

        if args.matrix_json and Path(args.matrix_json).exists():
            data = json.loads(Path(args.matrix_json).read_text(encoding="utf-8"))
        elif args.mock:
            data = load_mock_fixture("MA0139.1")
        else:
            raise ValueError("Must provide an existing --matrix-json file or pass --mock.")

        analysis = build_matrix_analysis(
            data,
            pseudocount=args.pseudocount,
            gc_content=args.gc_content,
        )

        if args.output_json:
            out_json = Path(args.output_json)
            out_json.parent.mkdir(parents=True, exist_ok=True)
            out_json.write_text(json.dumps(analysis, indent=2), encoding="utf-8")

        if args.output_html:
            html_doc = render_matrix_html_widget(analysis)
            out_html = Path(args.output_html)
            out_html.parent.mkdir(parents=True, exist_ok=True)
            out_html.write_text(html_doc, encoding="utf-8")

        if not args.output_json and not args.output_html:
            sys.stdout.write(json.dumps(analysis, indent=2) + "\n")
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
