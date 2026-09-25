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
"""End-to-End JASPAR Visualization and Sequence Scanning Pipeline Orchestrator.

Executes the complete workflow from a single command:
1. Resolves transcription factor gene symbol or JASPAR Matrix ID
2. Fetches PFM counts and TFClass structural metadata
3. Renders vector SVG and interactive HTML5 sequence logos
4. Computes PFM/PPM/PWM heatmaps and format exports (MEME, TRANSFAC, JASPAR, YAML)
5. Scans target FASTA or DNA sequence strings across forward (+) and reverse (-) strands
6. Generates the standalone interactive HTML5 dashboard (dashboard.html)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

from generate_dashboard import compile_dashboard_html
from jaspar_api import (
    ACADEMIC_CITATION_NOTICE,
    SCIENTIFIC_DISCLAIMER,
    PoliteHttpClient,
    JASPAR_API_BASE,
    load_mock_fixture,
    verify_workspace_license,
)
from visualize_logo import compute_logo_metrics, render_logo_html_widget, render_logo_svg
from visualize_matrix import build_matrix_analysis, render_matrix_html_widget
from visualize_scan import parse_fasta_or_raw, render_scan_html_widget, scan_dna_sequence


def execute_pipeline(
    name: Optional[str],
    matrix_id: Optional[str],
    tax_id: int,
    sequence: Optional[str],
    fasta: Optional[str],
    threshold: float,
    output_dir: str,
    mock: bool,
    motif_json: Optional[str] = None,
) -> Dict[str, Any]:
    """Runs the end-to-end analytical pipeline and writes all artifacts to output_dir."""
    verify_workspace_license(quiet=True)
    out_root = Path(output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    target_id = matrix_id or name or "MA0139.1"

    # Step 1 & 2: Resolve TF & Fetch Motif Profile (supports custom --motif-json, live JASPAR API, or --mock)
    if motif_json and Path(motif_json).exists():
        motif_data = json.loads(Path(motif_json).read_text(encoding="utf-8"))
    elif mock:
        motif_data = load_mock_fixture(target_id)
    else:
        client = PoliteHttpClient(qps=10.0)
        if matrix_id:
            motif_data = client.get_json(f"{JASPAR_API_BASE}/matrix/{matrix_id.strip()}/")
        else:
            tf_name = (name or "CTCF").strip()
            search_resp = client.get_json(f"{JASPAR_API_BASE}/matrix/?name={tf_name}&tax_id={tax_id}&version=latest")
            results = search_resp.get("results", [])
            if not results:
                raise ValueError(f"No JASPAR matrices found for name='{tf_name}', tax_id={tax_id}.")
            resolved_id = results[0]["matrix_id"]
            motif_data = client.get_json(f"{JASPAR_API_BASE}/matrix/{resolved_id}/")

    resolved_matrix_id = motif_data.get("matrix_id", "MA0139.1")
    resolved_name = motif_data.get("name", name or "CTCF")

    motif_json_path = out_root / "motif.json"
    motif_json_path.write_text(json.dumps(motif_data, indent=2), encoding="utf-8")

    # Step 3: Sequence Logo Generation
    logo_metrics = compute_logo_metrics(motif_data["pfm"], pseudocount=0.8, reverse_complement=False)
    logo_svg = render_logo_svg(logo_metrics, title=f"{resolved_matrix_id} ({resolved_name}) Sequence Logo", mode="bits")
    logo_html = render_logo_html_widget(motif_data, pseudocount=0.8, initial_revcomp=False)

    logo_svg_path = out_root / "logo.svg"
    logo_html_path = out_root / "logo.html"
    logo_svg_path.write_text(logo_svg, encoding="utf-8")
    logo_html_path.write_text(logo_html, encoding="utf-8")

    # Step 4: Matrix Heatmap Explorer
    matrix_analysis = build_matrix_analysis(motif_data, pseudocount=0.8, gc_content=0.50)
    matrix_json_path = out_root / "matrix.json"
    matrix_html_path = out_root / "matrix.html"
    matrix_json_path.write_text(json.dumps(matrix_analysis, indent=2), encoding="utf-8")
    matrix_html_path.write_text(render_matrix_html_widget(matrix_analysis), encoding="utf-8")

    # Step 5: DNA Sequence Scanning (processes user-provided --sequence or --fasta first)
    seq_name, seq_str = parse_fasta_or_raw(fasta, sequence, mock=True if (mock and not fasta and not sequence) else False)
    scan_results = scan_dna_sequence(
        matrix_data=motif_data,
        sequence=seq_str,
        sequence_name=seq_name,
        threshold=threshold,
    )
    scan_json_path = out_root / "scan_hits.json"
    scan_html_path = out_root / "scan_hits.html"
    scan_json_path.write_text(json.dumps(scan_results, indent=2), encoding="utf-8")
    scan_html_path.write_text(render_scan_html_widget(scan_results), encoding="utf-8")

    # Step 6: Unified Standalone Dashboard (dashboard.html) reflecting the processed user data
    dashboard_path = out_root / "dashboard.html"
    dashboard_html = compile_dashboard_html(
        motif_data=motif_data,
        metadata_data=motif_data,
        scan_data=scan_results,
        fasta_path=fasta,
        raw_sequence=sequence,
    )
    dashboard_path.write_text(dashboard_html, encoding="utf-8")

    summary = {
        "status": "success",
        "matrix_id": resolved_matrix_id,
        "name": resolved_name,
        "tax_id": tax_id,
        "motif_length": logo_metrics["length"],
        "consensus_forward": logo_metrics["consensus"],
        "total_information_bits": logo_metrics["total_bits"],
        "scanned_sequence": seq_name,
        "scanned_length_bp": scan_results["sequence_length"],
        "threshold": threshold,
        "total_binding_hits": scan_results["total_hits"],
        "artifacts": {
            "motif_json": str(motif_json_path),
            "logo_svg": str(logo_svg_path),
            "logo_html": str(logo_html_path),
            "matrix_json": str(matrix_json_path),
            "matrix_html": str(matrix_html_path),
            "scan_json": str(scan_json_path),
            "scan_html": str(scan_html_path),
            "dashboard_html": str(dashboard_path),
        },
        "citations": ACADEMIC_CITATION_NOTICE,
        "disclaimer": SCIENTIFIC_DISCLAIMER,
    }
    summary_path = out_root / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    sys.stdout.write(json.dumps(summary, indent=2) + "\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    """Constructs CLI argument parser for run_pipeline.py."""
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="End-to-end JASPAR TF resolution, logo rendering, sequence scanning, and dashboard pipeline.",
    )
    parser.add_argument("--name", help="Transcription factor name (e.g., CTCF, JUN, TP53, SOX2, MYC, FOXA1).")
    parser.add_argument("--matrix-id", help="Direct JASPAR Matrix ID (e.g., MA0139.1).")
    parser.add_argument("--motif-json", help="Optional path to a custom user PFM motif JSON file.")
    parser.add_argument("--tax-id", type=int, default=9606, help="NCBI Taxonomy ID (default: 9606).")
    parser.add_argument("--sequence", help="User-provided DNA sequence string (ACGT...) to scan.")
    parser.add_argument("--fasta", help="User-provided FASTA file path to scan.")
    parser.add_argument("--threshold", type=float, default=0.80, help="Relative score scan threshold (default: 0.80).")
    parser.add_argument("--output-dir", required=True, help="Target directory for generated artifacts and dashboard.html.")
    parser.add_argument("--mock", action="store_true", help="Execute offline using mock motif fixtures.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for run_pipeline.py."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        execute_pipeline(
            name=args.name,
            matrix_id=args.matrix_id,
            tax_id=args.tax_id,
            sequence=args.sequence,
            fasta=args.fasta,
            threshold=args.threshold,
            output_dir=args.output_dir,
            mock=args.mock,
            motif_json=args.motif_json,
        )
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
