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
"""Information-Theoretic Sequence Logo Engine (SVG and Interactive HTML5).

Computes Position Probability Matrices (PPM), Shannon entropy (H_i), Schneider
small-sample error correction (e_n), total information content (R_i in 0..2 bits),
base letter heights, and IUPAC consensus nucleotide symbols. Renders vector
SVG <path> glyphs without external font dependencies.
"""

from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

from jaspar_api import (
    ACADEMIC_CITATION_NOTICE,
    SCIENTIFIC_DISCLAIMER,
    load_mock_fixture,
    validate_pfm,
    verify_workspace_license,
)

BASE_COLORS: Dict[str, str] = {
    "A": "#10B981",  # Green: rgb(16, 185, 129)
    "C": "#3B82F6",  # Blue:  rgb(59, 130, 246)
    "G": "#F59E0B",  # Orange: rgb(245, 158, 11)
    "T": "#EF4444",  # Red:   rgb(239, 68, 68)
}

# Vector SVG <path> glyph definitions normalized to a 100 x 100 bounding box (0..100 X, 0..100 Y)
BASE_SVG_PATHS: Dict[str, str] = {
    "A": (
        "M 50 4 L 94 96 L 73 96 L 63 72 L 37 72 L 27 96 L 6 96 Z "
        "M 50 28 L 42 56 L 58 56 Z"
    ),
    "C": (
        "M 88 26 C 78 10, 62 4, 46 4 C 20 4, 6 24, 6 50 "
        "C 6 76, 20 96, 46 96 C 64 96, 79 88, 88 72 L 70 63 "
        "C 64 74, 56 79, 46 79 C 31 79, 24 67, 24 50 "
        "C 24 33, 31 21, 46 21 C 56 21, 64 26, 70 36 Z"
    ),
    "G": (
        "M 90 26 C 80 10, 64 4, 47 4 C 20 4, 6 24, 6 50 "
        "C 6 76, 20 96, 47 96 C 68 96, 86 84, 92 64 L 92 46 "
        "L 50 46 L 50 62 L 73 62 C 70 73, 60 79, 47 79 "
        "C 31 79, 24 67, 24 50 C 24 33, 31 21, 47 21 "
        "C 58 21, 67 26, 72 36 Z"
    ),
    "T": (
        "M 8 4 L 92 4 L 92 22 L 60 22 L 60 96 L 40 96 L 40 22 L 8 22 Z"
    ),
}

IUPAC_PAIR_MAP: Dict[Tuple[str, str], str] = {
    ("A", "G"): "R",
    ("C", "T"): "Y",
    ("C", "G"): "S",
    ("A", "T"): "W",
    ("G", "T"): "K",
    ("A", "C"): "M",
}

IUPAC_TRIPLE_MAP: Dict[Tuple[str, str, str], str] = {
    ("C", "G", "T"): "B",
    ("A", "G", "T"): "D",
    ("A", "C", "T"): "H",
    ("A", "C", "G"): "V",
}


def reverse_complement_pfm(pfm: Dict[str, List[float]]) -> Dict[str, List[float]]:
    """Computes the reverse complement of a PFM matrix (reversing columns and swapping A<->T, C<->G)."""
    validate_pfm(pfm)
    return {
        "A": list(reversed(pfm["T"])),
        "C": list(reversed(pfm["G"])),
        "G": list(reversed(pfm["C"])),
        "T": list(reversed(pfm["A"])),
    }


def derive_iupac_symbol(probs: Dict[str, float]) -> str:
    """Derives the IUPAC consensus nucleotide symbol for a single column per Section 5.3.6."""
    sorted_bases = sorted(probs.items(), key=lambda kv: kv[1], reverse=True)
    b1, p1 = sorted_bases[0]
    b2, p2 = sorted_bases[1]
    b3, p3 = sorted_bases[2]

    # Single base > 50% and >= twice the second most frequent base
    if p1 > 0.50 and p1 >= 2.0 * p2:
        return b1
    # Two bases sum > 75%
    if (p1 + p2) > 0.75:
        pair_key = tuple(sorted((b1, b2)))
        return IUPAC_PAIR_MAP.get(pair_key, "N")  # type: ignore
    # Three bases sum > 85%
    if (p1 + p2 + p3) > 0.85:
        triple_key = tuple(sorted((b1, b2, b3)))
        return IUPAC_TRIPLE_MAP.get(triple_key, "N")  # type: ignore
    return "N"


def compute_logo_metrics(
    pfm: Dict[str, List[float]],
    pseudocount: float = 0.8,
    reverse_complement: bool = False,
) -> Dict[str, Any]:
    """Computes Shannon entropy H_i, small-sample error e_n, R_i, base heights, and IUPAC strip."""
    if pseudocount <= 0.0 or pseudocount > 100.0:
        raise ValueError(f"Pseudocount must be in (0.0, 100.0], received {pseudocount}.")

    working_pfm = reverse_complement_pfm(pfm) if reverse_complement else pfm
    length = validate_pfm(working_pfm)

    columns: List[Dict[str, Any]] = []
    consensus_chars: List[str] = []
    total_bits = 0.0

    for i in range(length):
        n_i = sum(float(working_pfm[b][i]) for b in ("A", "C", "G", "T"))
        denom = n_i + 4.0 * pseudocount
        probs = {b: (float(working_pfm[b][i]) + pseudocount) / denom for b in ("A", "C", "G", "T")}

        # 1. Shannon Entropy H_i = - sum(p * log2(p))
        h_i = 0.0
        for b in ("A", "C", "G", "T"):
            p = probs[b]
            if p > 0.0:
                h_i -= p * math.log2(p)

        # 2. Small-Sample Error Correction e_n = 3 / (2 * ln(2) * N_i)
        if n_i > 0.0:
            e_n = 3.0 / (2.0 * math.log(2.0) * n_i)
        else:
            e_n = 0.0

        # 3. Total Information Content R_i = max(0, 2 - (H_i + e_n))
        r_i = max(0.0, 2.0 - (h_i + e_n))
        total_bits += r_i

        # 4. Letter Heights = PPM[b, i] * R_i
        heights = {b: probs[b] * r_i for b in ("A", "C", "G", "T")}

        # 5. Stack order: ascending by height (smallest on bottom, largest on top)
        stack = sorted(
            [
                {
                    "base": b,
                    "prob": round(probs[b], 6),
                    "bits": round(heights[b], 6),
                }
                for b in ("A", "C", "G", "T")
            ],
            key=lambda item: (item["bits"], item["prob"], item["base"]),
        )

        iupac = derive_iupac_symbol(probs)
        consensus_chars.append(iupac)

        columns.append(
            {
                "position": i + 1,
                "n_i": round(n_i, 2),
                "entropy_h_i": round(h_i, 6),
                "error_e_n": round(e_n, 6),
                "info_r_i": round(r_i, 6),
                "probs": {b: round(probs[b], 6) for b in ("A", "C", "G", "T")},
                "heights": {b: round(heights[b], 6) for b in ("A", "C", "G", "T")},
                "stack": stack,
                "iupac": iupac,
            }
        )

    return {
        "length": length,
        "pseudocount": pseudocount,
        "reverse_complement": reverse_complement,
        "total_bits": round(total_bits, 4),
        "consensus": "".join(consensus_chars),
        "columns": columns,
    }


def render_logo_svg(
    metrics: Dict[str, Any],
    title: str = "JASPAR Sequence Logo",
    mode: str = "bits",
) -> str:
    """Renders a standalone vector SVG sequence logo using pure <path> glyphs."""
    length = metrics["length"]
    col_width = 42.0
    left_margin = 68.0
    right_margin = 28.0
    top_margin = 52.0
    plot_height = 220.0
    bottom_margin = 76.0

    svg_width = int(left_margin + length * col_width + right_margin)
    svg_height = int(top_margin + plot_height + bottom_margin)
    max_y_val = 2.0 if mode == "bits" else 1.0
    y_label = "Information Content (bits)" if mode == "bits" else "Position Probability"
    y_ticks = [0.0, 0.5, 1.0, 1.5, 2.0] if mode == "bits" else [0.0, 0.25, 0.5, 0.75, 1.0]

    orientation_tag = "Reverse Complement (3' <- 5')" if metrics["reverse_complement"] else "Forward (5' -> 3')"
    escaped_title = html.escape(f"{title} [{orientation_tag}]")

    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_width} {svg_height}" '
        f'width="100%" height="{svg_height}" role="img" aria-label="{escaped_title}" '
        'style="background:#FFFFFF; font-family: ui-sans-serif, system-ui, -apple-system, sans-serif;">',
        f'  <title>{escaped_title}</title>',
        f'  <rect x="0" y="0" width="{svg_width}" height="{svg_height}" fill="#FFFFFF" rx="8"/>',
        f'  <text x="{left_margin}" y="28" font-size="15" font-weight="700" fill="#0F172A">{escaped_title}</text>',
        f'  <text x="{svg_width - right_margin}" y="28" font-size="12" font-weight="600" fill="#475569" text-anchor="end">'
        f'Total Info: {metrics["total_bits"]:.2f} bits | Consensus: {metrics["consensus"]}</text>',
    ]

    # Y-axis grid lines and labels
    for tick in y_ticks:
        frac = tick / max_y_val
        y_pos = top_margin + plot_height * (1.0 - frac)
        stroke_color = "#CBD5E1" if tick == 0.0 else "#E2E8F0"
        dash = "" if tick in (0.0, max_y_val) else ' stroke-dasharray="4,4"'
        parts.append(
            f'  <line x1="{left_margin}" y1="{y_pos:.2f}" x2="{svg_width - right_margin}" y2="{y_pos:.2f}" '
            f'stroke="{stroke_color}" stroke-width="1"{dash}/>'
        )
        parts.append(
            f'  <text x="{left_margin - 10}" y="{y_pos + 4:.2f}" font-size="11" fill="#475569" text-anchor="end">'
            f'{tick:.1f}</text>'
        )

    # Y-axis spine and rotated title
    parts.append(
        f'  <line x1="{left_margin}" y1="{top_margin}" x2="{left_margin}" y2="{top_margin + plot_height}" '
        'stroke="#64748B" stroke-width="1.5"/>'
    )
    mid_y = top_margin + plot_height / 2.0
    parts.append(
        f'  <text x="18" y="{mid_y:.2f}" font-size="12" font-weight="600" fill="#334155" '
        f'text-anchor="middle" transform="rotate(-90, 18, {mid_y:.2f})">{y_label}</text>'
    )

    # Render stacked letter columns
    baseline_y = top_margin + plot_height
    for col_idx, col in enumerate(metrics["columns"]):
        x_left = left_margin + col_idx * col_width
        x_center = x_left + col_width / 2.0
        current_y = baseline_y

        for item in col["stack"]:
            base = item["base"]
            val = item["bits"] if mode == "bits" else item["prob"]
            if val <= 0.002:
                continue
            glyph_h = (val / max_y_val) * plot_height
            top_y = current_y - glyph_h
            scale_x = (col_width - 4.0) / 100.0
            scale_y = glyph_h / 100.0
            path_d = BASE_SVG_PATHS[base]
            color = BASE_COLORS[base]
            tooltip = (
                f"Pos {col['position']} ({base}): {item['bits']:.3f} bits | "
                f"P={item['prob']:.3f} | R_i={col['info_r_i']:.3f} bits"
            )
            parts.append(
                f'  <g transform="translate({x_left + 2.0:.2f}, {top_y:.2f}) scale({scale_x:.4f}, {scale_y:.4f})">'
                f'<title>{tooltip}</title>'
                f'<path d="{path_d}" fill="{color}"/>'
                f'</g>'
            )
            current_y = top_y

        # Position index tick label
        parts.append(
            f'  <text x="{x_center:.2f}" y="{baseline_y + 18:.2f}" font-size="11" font-weight="600" '
            f'fill="#334155" text-anchor="middle">{col["position"]}</text>'
        )

        # Bottom IUPAC Consensus Strip Badge
        iupac_char = col["iupac"]
        badge_bg = BASE_COLORS.get(iupac_char, "#64748B")
        parts.append(
            f'  <rect x="{x_left + 3:.2f}" y="{baseline_y + 28:.2f}" width="{col_width - 6:.2f}" height="24" '
            f'rx="4" fill="{badge_bg}" fill-opacity="0.15" stroke="{badge_bg}" stroke-width="1"/>'
        )
        parts.append(
            f'  <text x="{x_center:.2f}" y="{baseline_y + 44:.2f}" font-size="12" font-weight="700" '
            f'fill="#0F172A" text-anchor="middle">{iupac_char}</text>'
        )

    parts.append(
        f'  <text x="{left_margin}" y="{svg_height - 8}" font-size="11" font-weight="600" fill="#64748B">'
        'IUPAC Consensus Strip (Position 1 -> L)</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def render_logo_html_widget(
    data: Dict[str, Any],
    pseudocount: float = 0.8,
    initial_revcomp: bool = False,
    title: Optional[str] = None,
) -> str:
    """Generates a standalone interactive HTML5 Sequence Logo widget with client-side toggles."""
    matrix_id = data.get("matrix_id", "MA0139.1")
    name = data.get("name", "TF")
    header_title = title or f"{matrix_id} ({name}) Sequence Logo"

    fwd_metrics = compute_logo_metrics(data["pfm"], pseudocount=pseudocount, reverse_complement=False)
    rev_metrics = compute_logo_metrics(data["pfm"], pseudocount=pseudocount, reverse_complement=True)

    svg_fwd_bits = render_logo_svg(fwd_metrics, title=header_title, mode="bits")
    svg_rev_bits = render_logo_svg(rev_metrics, title=header_title, mode="bits")
    svg_fwd_prob = render_logo_svg(fwd_metrics, title=header_title, mode="prob")
    svg_rev_prob = render_logo_svg(rev_metrics, title=header_title, mode="prob")

    state_payload = {
        "matrix_id": matrix_id,
        "name": name,
        "initial_revcomp": initial_revcomp,
        "svg": {
            "fwd_bits": svg_fwd_bits,
            "rev_bits": svg_rev_bits,
            "fwd_prob": svg_fwd_prob,
            "rev_prob": svg_rev_prob,
        },
        "consensus": {
            "fwd": fwd_metrics["consensus"],
            "rev": rev_metrics["consensus"],
        },
    }

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>{html.escape(header_title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, sans-serif; background: #F8FAFC; color: #0F172A; margin: 24px; }}
    .card {{ background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 20px; max-width: 1040px; margin: 0 auto; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
    .toolbar {{ display: flex; gap: 10px; flex-wrap: wrap; align-items: center; justify-content: space-between; margin-bottom: 16px; }}
    .btn {{ background: #0F172A; color: #FFFFFF; border: none; border-radius: 6px; padding: 8px 14px; font-size: 13px; font-weight: 600; cursor: pointer; }}
    .btn.secondary {{ background: #E2E8F0; color: #0F172A; }}
    .consensus-box {{ background: #F1F5F9; padding: 10px 14px; border-radius: 6px; font-family: monospace; font-size: 14px; margin-top: 12px; }}
    .footer {{ font-size: 12px; color: #475569; margin-top: 16px; border-top: 1px solid #E2E8F0; padding-top: 12px; line-height: 1.5; }}
  </style>
</head>
<body>
  <div class="card">
    <div class="toolbar">
      <div>
        <strong>{html.escape(header_title)}</strong>
      </div>
      <div>
        <button class="btn" id="toggleStrandBtn" onclick="toggleStrand()">Toggle Reverse Complement</button>
        <button class="btn secondary" id="toggleModeBtn" onclick="toggleMode()">Switch Y-Axis (Bits / Prob)</button>
      </div>
    </div>
    <div id="logoContainer"></div>
    <div class="consensus-box">
      IUPAC Consensus: <strong id="consensusText"></strong>
    </div>
    <div class="footer">
      <p><strong>Academic Citation:</strong> Ovek Baydar, D. et al. JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models. Nucleic Acids Research 54, D184-D193 (2026). doi:10.1093/nar/gkaf1209. Licensed under CC BY 4.0.</p>
      <p><strong>Scientific Disclaimer:</strong> {html.escape(SCIENTIFIC_DISCLAIMER)}</p>
    </div>
  </div>
  <script>
    const LOGO_STATE = {json.dumps(state_payload)};
    let isRevComp = Boolean(LOGO_STATE.initial_revcomp);
    let isProbMode = false;

    function renderLogo() {{
      const key = (isRevComp ? "rev_" : "fwd_") + (isProbMode ? "prob" : "bits");
      document.getElementById("logoContainer").innerHTML = LOGO_STATE.svg[key];
      document.getElementById("consensusText").textContent = isRevComp ? LOGO_STATE.consensus.rev : LOGO_STATE.consensus.fwd;
    }}
    function toggleStrand() {{ isRevComp = !isRevComp; renderLogo(); }}
    function toggleMode() {{ isProbMode = !isProbMode; renderLogo(); }}
    renderLogo();
  </script>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    """Constructs CLI argument parser for visualize_logo.py."""
    parser = argparse.ArgumentParser(
        prog="visualize_logo.py",
        description="Compute Shannon information content and generate vector SVG / HTML5 sequence logos.",
    )
    parser.add_argument("--matrix-json", help="Path to JSON file containing matrix PFM data.")
    parser.add_argument("--output-svg", help="Path to save standalone .svg sequence logo.")
    parser.add_argument("--output-html", help="Path to save interactive HTML5 sequence logo widget.")
    parser.add_argument("--pseudocount", type=float, default=0.8, help="Per-base pseudocount (default: 0.8).")
    parser.add_argument("--title", help="Title string for the logo graphic.")
    parser.add_argument("--reverse-complement", action="store_true", help="Generate the reverse complement logo.")
    parser.add_argument("--mock", action="store_true", help="Use embedded mock motif if --matrix-json is omitted.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for visualize_logo.py."""
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

        title = args.title or f"{data.get('matrix_id', 'MA0139.1')} {data.get('name', 'CTCF')} Sequence Logo"
        metrics = compute_logo_metrics(
            data["pfm"],
            pseudocount=args.pseudocount,
            reverse_complement=args.reverse_complement,
        )

        if args.output_svg:
            svg_text = render_logo_svg(metrics, title=title, mode="bits")
            out_svg = Path(args.output_svg)
            out_svg.parent.mkdir(parents=True, exist_ok=True)
            out_svg.write_text(svg_text, encoding="utf-8")

        if args.output_html:
            html_text = render_logo_html_widget(
                data,
                pseudocount=args.pseudocount,
                initial_revcomp=args.reverse_complement,
                title=title,
            )
            out_html = Path(args.output_html)
            out_html.parent.mkdir(parents=True, exist_ok=True)
            out_html.write_text(html_text, encoding="utf-8")

        if not args.output_svg and not args.output_html:
            sys.stdout.write(json.dumps(metrics, indent=2) + "\n")
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
