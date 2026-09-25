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
"""Sliding-Window DNA Sequence Scanner and Dual-Strand Motif Matcher.

Scans FASTA files or nucleotide strings (up to 100kb) across forward (+) and
reverse-complement (-) strands using log-odds Position Weight Matrices (PWMs).
Computes raw scores (bits), normalized relative scores (S_rel in [0.0, 1.0]),
and categorical significance tiers (Optimal, High, Moderate).
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

from jaspar_api import (
    ACADEMIC_CITATION_NOTICE,
    SCIENTIFIC_DISCLAIMER,
    compute_ppm_and_pwm,
    find_skill_root,
    load_mock_fixture,
    verify_workspace_license,
)

MAX_SEQUENCE_LENGTH = 100000
COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement_seq(seq: str) -> str:
    """Returns the 5'->3' reverse complement of a nucleotide string."""
    return seq.translate(COMPLEMENT_TABLE)[::-1].upper()


def parse_fasta_or_raw(fasta_path: Optional[str], raw_seq: Optional[str], mock: bool = False) -> Tuple[str, str]:
    """Loads sequence header and clean nucleotide string from FASTA path, raw sequence, or mock fixture."""
    if raw_seq:
        cleaned = "".join(raw_seq.split()).upper()
        return ("USER_PROVIDED_SEQUENCE", cleaned)

    target_path: Optional[Path] = None
    if fasta_path and Path(fasta_path).exists():
        target_path = Path(fasta_path)
    elif mock:
        target_path = find_skill_root() / "sample_data" / "sample_enhancer.fasta"

    if not target_path or not target_path.exists():
        raise ValueError("Must provide an existing --fasta file, --sequence string, or pass --mock.")

    lines = target_path.read_text(encoding="utf-8").splitlines()
    header = "SEQUENCE"
    seq_parts: List[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if not seq_parts:
                header = stripped[1:].split()[0]
            else:
                break
        else:
            seq_parts.append(stripped.upper())

    full_seq = "".join(seq_parts)
    return (header, full_seq)


def classify_significance(rel_score: float) -> str:
    """Assigns categorical rank based on relative affinity score."""
    if rel_score >= 0.95:
        return "Optimal"
    if rel_score >= 0.85:
        return "High"
    if rel_score >= 0.80:
        return "Moderate"
    return "Candidate"


def scan_dna_sequence(
    matrix_data: Dict[str, Any],
    sequence: str,
    sequence_name: str = "SEQUENCE",
    threshold: float = 0.80,
    pseudocount: float = 0.8,
    gc_content: float = 0.50,
) -> Dict[str, Any]:
    """Scans sequence across forward (+) and reverse (-) strands using sliding-window PWM log-odds."""
    if not (0.0 <= threshold <= 1.0):
        raise ValueError(f"Threshold must be within [0.0, 1.0], received {threshold}.")
    seq_clean = "".join(sequence.split()).upper()
    seq_len = len(seq_clean)
    if seq_len > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence length ({seq_len} bp) exceeds the 100,000 bp (100kb) safety limit. "
            "Chunk genomic sequences into <= 100kb windows before scanning."
        )
    invalid_chars = set(seq_clean) - {"A", "C", "G", "T", "N"}
    if invalid_chars:
        raise ValueError(f"Sequence contains unsupported non-DNA characters: {sorted(invalid_chars)}.")

    calc = compute_ppm_and_pwm(
        matrix_data["pfm"],
        pseudocount=pseudocount,
        gc_content=gc_content,
    )
    pwm = calc["pwm"]
    motif_len = calc["length"]

    # Compute theoretical extreme scores S_min and S_max
    s_max = sum(max(pwm[b][i] for b in ("A", "C", "G", "T")) for i in range(motif_len))
    s_min = sum(min(pwm[b][i] for b in ("A", "C", "G", "T")) for i in range(motif_len))
    score_span = max(s_max - s_min, 1e-9)

    hits: List[Dict[str, Any]] = []
    warnings: List[str] = []

    if seq_len < motif_len:
        warnings.append(
            f"Input sequence length ({seq_len} bp) is shorter than motif length ({motif_len} bp); 0 windows scanned."
        )
    else:
        for offset in range(seq_len - motif_len + 1):
            fwd_kmer = seq_clean[offset : offset + motif_len]
            rev_kmer = reverse_complement_seq(fwd_kmer)

            for strand, kmer in (("+", fwd_kmer), ("-", rev_kmer)):
                raw_score = 0.0
                for pos_idx, base in enumerate(kmer):
                    if base in ("A", "C", "G", "T"):
                        raw_score += pwm[base][pos_idx]
                    else:
                        # Neutral log-odds weight (0.0 bits) for ambiguous base 'N'
                        raw_score += 0.0

                rel_score = (raw_score - s_min) / score_span
                rel_score = max(0.0, min(1.0, rel_score))

                if rel_score >= threshold:
                    hits.append(
                        {
                            "start": offset + 1,
                            "end": offset + motif_len,
                            "strand": strand,
                            "matched_sequence": kmer,
                            "raw_score": round(raw_score, 4),
                            "relative_score": round(rel_score, 4),
                            "significance": classify_significance(rel_score),
                        }
                    )

    # Sort hits by relative_score descending, then by genomic coordinate ascending
    hits.sort(key=lambda h: (-h["relative_score"], h["start"], h["strand"]))

    return {
        "matrix_id": matrix_data.get("matrix_id", "MA0139.1"),
        "name": matrix_data.get("name", "CTCF"),
        "motif_length": motif_len,
        "sequence_name": sequence_name,
        "sequence_length": seq_len,
        "sequence": seq_clean,
        "threshold": threshold,
        "pseudocount": pseudocount,
        "gc_content": gc_content,
        "s_min": round(s_min, 4),
        "s_max": round(s_max, 4),
        "total_hits": len(hits),
        "warnings": warnings,
        "hits": hits,
        "citations": ACADEMIC_CITATION_NOTICE,
        "disclaimer": SCIENTIFIC_DISCLAIMER,
    }


def render_locus_map_svg(scan_result: Dict[str, Any]) -> str:
    """Renders a vector SVG genomic locus distribution track displaying forward and reverse hits."""
    seq_len = max(1, scan_result["sequence_length"])
    width = 920
    height = 170
    left_pad = 70
    right_pad = 30
    track_w = width - left_pad - right_pad
    mid_y = 88

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="{height}" '
        'style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; font-family: system-ui, sans-serif;">',
        f'  <text x="{left_pad}" y="24" font-size="13" font-weight="700" fill="#0F172A">'
        f'Binding Site Locus Map: {html.escape(scan_result["sequence_name"])} (1 .. {seq_len} bp)</text>',
        f'  <line x1="{left_pad}" y1="{mid_y}" x2="{width - right_pad}" y2="{mid_y}" stroke="#334155" stroke-width="2"/>',
        f'  <text x="{left_pad - 10}" y="{mid_y - 26}" font-size="11" font-weight="600" fill="#2563EB" text-anchor="end">+ Strand</text>',
        f'  <text x="{left_pad - 10}" y="{mid_y + 34}" font-size="11" font-weight="600" fill="#DC2626" text-anchor="end">- Strand</text>',
    ]

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = left_pad + frac * track_w
        bp_coord = int(round(1 + frac * (seq_len - 1)))
        parts.append(f'  <line x1="{x:.1f}" y1="{mid_y - 6}" x2="{x:.1f}" y2="{mid_y + 6}" stroke="#475569" stroke-width="1"/>')
        parts.append(f'  <text x="{x:.1f}" y="{height - 14}" font-size="10" fill="#64748B" text-anchor="middle">{bp_coord} bp</text>')

    for hit in scan_result["hits"]:
        x = left_pad + ((hit["start"] - 1) / max(seq_len, 1)) * track_w
        is_fwd = hit["strand"] == "+"
        color = "#2563EB" if is_fwd else "#DC2626"
        stem_h = 18 + (hit["relative_score"] - 0.5) * 44
        y_end = (mid_y - stem_h) if is_fwd else (mid_y + stem_h)
        parts.append(f'  <line x1="{x:.1f}" y1="{mid_y}" x2="{x:.1f}" y2="{y_end:.1f}" stroke="{color}" stroke-width="2.5"/>')
        parts.append(
            f'  <circle cx="{x:.1f}" cy="{y_end:.1f}" r="5.5" fill="{color}">'
            f'<title>{hit["strand"]} Strand [{hit["start"]}..{hit["end"]}]: {hit["matched_sequence"]} '
            f'(Rel={hit["relative_score"]:.3f}, Raw={hit["raw_score"]:.2f} bits, {hit["significance"]})</title>'
            f'</circle>'
        )

    parts.append("</svg>")
    return "\n".join(parts)


def render_scan_html_widget(scan_result: Dict[str, Any]) -> str:
    """Renders an interactive HTML5 sequence scanner report with locus map and filterable hit table."""
    title = f"{scan_result['matrix_id']} ({scan_result['name']}) Sequence Scan Report"
    locus_svg = render_locus_map_svg(scan_result)

    rows_html = []
    for idx, h in enumerate(scan_result["hits"], start=1):
        rows_html.append(
            f"<tr>"
            f"<td>{idx}</td>"
            f"<td>{h['start']}</td>"
            f"<td>{h['end']}</td>"
            f"<td><strong>{h['strand']}</strong></td>"
            f"<td style='font-family:monospace'>{html.escape(h['matched_sequence'])}</td>"
            f"<td>{h['raw_score']:.3f}</td>"
            f"<td>{(h['relative_score'] * 100):.1f}% ({h['relative_score']:.3f})</td>"
            f"<td>{html.escape(h['significance'])}</td>"
            f"</tr>"
        )
    tbody = "\n".join(rows_html) if rows_html else "<tr><td colspan='8'>No binding sites above threshold.</td></tr>"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, sans-serif; background: #F8FAFC; color: #0F172A; margin: 24px; }}
    .card {{ background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 20px; max-width: 1100px; margin: 0 auto; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 16px; font-size: 13px; }}
    th, td {{ border: 1px solid #E2E8F0; padding: 8px 10px; text-align: left; }}
    th {{ background: #F1F5F9; font-weight: 700; }}
    .footer {{ font-size: 12px; color: #475569; margin-top: 16px; border-top: 1px solid #E2E8F0; padding-top: 12px; }}
  </style>
</head>
<body>
  <div class="card">
    <h2>{html.escape(title)}</h2>
    <p>Sequence: <strong>{html.escape(scan_result['sequence_name'])}</strong> ({scan_result['sequence_length']} bp) |
       Threshold: <strong>{scan_result['threshold']:.2f}</strong> |
       Total Hits: <strong>{scan_result['total_hits']}</strong></p>
    {locus_svg}
    <table>
      <thead>
        <tr>
          <th>#</th><th>Start</th><th>End</th><th>Strand</th><th>Matched Sequence (5'->3')</th>
          <th>Raw Score (bits)</th><th>Relative Score</th><th>Significance</th>
        </tr>
      </thead>
      <tbody>
        {tbody}
      </tbody>
    </table>
    <div class="footer">
      <p><strong>Academic Citation:</strong> Ovek Baydar, D. et al. JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models. Nucleic Acids Research 54, D184-D193 (2026). doi:10.1093/nar/gkaf1209. Licensed under CC BY 4.0.</p>
      <p><strong>Scientific Disclaimer:</strong> {html.escape(SCIENTIFIC_DISCLAIMER)}</p>
    </div>
  </div>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    """Constructs CLI argument parser for visualize_scan.py."""
    parser = argparse.ArgumentParser(
        prog="visualize_scan.py",
        description="Scan DNA sequences or FASTA files on both strands using PWM log-odds scoring.",
    )
    parser.add_argument("--matrix-json", help="Path to matrix JSON file.")
    parser.add_argument("--fasta", help="Path to FASTA file to scan.")
    parser.add_argument("--sequence", help="Raw nucleotide sequence string (ACGT...) to scan.")
    parser.add_argument("--threshold", type=float, default=0.80, help="Relative score cutoff in [0.0, 1.0] (default: 0.80).")
    parser.add_argument("--pseudocount", type=float, default=0.8, help="Per-base pseudocount (default: 0.8).")
    parser.add_argument("--gc-content", type=float, default=0.50, help="Background GC frequency (default: 0.50).")
    parser.add_argument("--output-json", help="Path to save matched binding sites in structured JSON format.")
    parser.add_argument("--output-html", help="Path to save interactive match table and locus distribution map.")
    parser.add_argument("--mock", action="store_true", help="Offline mock scanning against sample_data/sample_enhancer.fasta.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for visualize_scan.py."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        verify_workspace_license(quiet=True)

        if args.matrix_json and Path(args.matrix_json).exists():
            matrix_data = json.loads(Path(args.matrix_json).read_text(encoding="utf-8"))
        elif args.mock:
            matrix_data = load_mock_fixture("MA0139.1")
        else:
            raise ValueError("Must provide an existing --matrix-json file or pass --mock.")

        seq_name, seq_str = parse_fasta_or_raw(args.fasta, args.sequence, mock=args.mock)
        scan_result = scan_dna_sequence(
            matrix_data=matrix_data,
            sequence=seq_str,
            sequence_name=seq_name,
            threshold=args.threshold,
            pseudocount=args.pseudocount,
            gc_content=args.gc_content,
        )

        if args.output_json:
            out_json = Path(args.output_json)
            out_json.parent.mkdir(parents=True, exist_ok=True)
            out_json.write_text(json.dumps(scan_result, indent=2), encoding="utf-8")

        if args.output_html:
            html_doc = render_scan_html_widget(scan_result)
            out_html = Path(args.output_html)
            out_html.parent.mkdir(parents=True, exist_ok=True)
            out_html.write_text(html_doc, encoding="utf-8")

        if not args.output_json and not args.output_html:
            sys.stdout.write(json.dumps(scan_result, indent=2) + "\n")
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
