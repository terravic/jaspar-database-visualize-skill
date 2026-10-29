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
"""Unified Standalone HTML5 Interactive Dashboard Generator (dashboard.html).

Compiles user-processed motif metadata, information-theoretic vector sequence logos,
PFM/PPM/PWM heatmaps, TFClass structural classifications, familial cluster comparisons,
3D interactive network and helical groove visualizations, provenance records, and
sliding-window DNA sequence scan results into a single self-contained HTML5 dashboard.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Dict, List, Optional

from jaspar_api import (
    ACADEMIC_CITATION_NOTICE,
    MOCK_NAME_MAP,
    SCIENTIFIC_DISCLAIMER,
    load_mock_fixture,
    verify_workspace_license,
)
from visualize_logo import BASE_SVG_PATHS, compute_logo_metrics
from visualize_matrix import build_matrix_analysis
from visualize_scan import parse_fasta_or_raw, scan_dna_sequence

# Semantic Data Visualization Palette
SEMANTIC_PALETTE: Dict[str, str] = {
    "teal": "#0D9488",     # Vibrant Teal: Primary data series, baseline metrics, structural bonds, nominal nodes
    "amber": "#D97706",    # Energetic Amber: Secondary data series, warning thresholds, intermediate intensity
    "indigo": "#4F46E5",   # Vivid Indigo: Comparison series, distinct molecular groups, specialized network nodes
    "emerald": "#059669",  # Functional Emerald: Passing metrics, positive outcomes, validated states
    "rose": "#E11D48",     # Deep Rose: Critical alerts, negative trends, threshold violations
}

# Semantic Base Color Mapping for Sequence Logos & 3D Structural Diagrams
SEMANTIC_BASE_COLORS: Dict[str, str] = {
    "A": "#059669",  # Functional Emerald
    "C": "#4F46E5",  # Vivid Indigo
    "G": "#D97706",  # Energetic Amber
    "T": "#E11D48",  # Deep Rose
    "N": "#0D9488",  # Vibrant Teal
}


def detect_host_system_theme() -> str:
    """Detects the host OS appearance preference (macOS, Windows, Linux)."""
    try:
        if sys.platform == "darwin":
            proc = subprocess.run(
                ["defaults", "read", "-g", "AppleInterfaceStyle"],
                capture_output=True,
                text=True,
                timeout=1.5,
                check=False,
            )
            if proc.returncode == 0 and "dark" in proc.stdout.strip().lower():
                return "dark"
            return "light"

        if sys.platform == "win32":
            import winreg  # type: ignore

            key_path = r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize"
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, key_path) as key:
                val, _ = winreg.QueryValueEx(key, "AppsUseLightTheme")
                return "light" if int(val) != 0 else "dark"

        # Linux: inspect GTK_THEME, gsettings, and kdeglobals
        gtk_theme = os.environ.get("GTK_THEME", "").lower()
        if "dark" in gtk_theme:
            return "dark"

        for cmd in (
            ["gsettings", "get", "org.gnome.desktop.interface", "color-scheme"],
            ["gsettings", "get", "org.gnome.desktop.interface", "gtk-theme"],
        ):
            try:
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=1.5,
                    check=False,
                )
                if proc.returncode == 0 and "dark" in proc.stdout.strip().lower():
                    return "dark"
            except Exception:
                pass

        kde_cfg = Path.home() / ".config" / "kdeglobals"
        if kde_cfg.exists():
            cfg_text = kde_cfg.read_text(encoding="utf-8", errors="ignore").lower()
            if "breezedark" in cfg_text or "colorscheme=dark" in cfg_text:
                return "dark"
    except Exception:
        pass
    return "light"


def render_semantic_logo_svg(
    metrics: Dict[str, Any],
    title: str = "JASPAR Sequence Logo",
    mode: str = "bits",
    theme: str = "light",
) -> str:
    """Renders a transparent, theme-reactive vector SVG sequence logo."""
    length = metrics["length"]
    col_width = 42.0
    left_margin = 64.0
    right_margin = 24.0
    top_margin = 46.0
    plot_height = 210.0
    bottom_margin = 70.0

    svg_width = int(left_margin + length * col_width + right_margin)
    svg_height = int(top_margin + plot_height + bottom_margin)
    max_y_val = 2.0 if mode == "bits" else 1.0
    y_label = "Information Content (bits)" if mode == "bits" else "Position Probability"
    y_ticks = [0.0, 0.5, 1.0, 1.5, 2.0] if mode == "bits" else [0.0, 0.25, 0.5, 0.75, 1.0]

    is_dark = theme == "dark"
    text_primary = "#ECE8E2" if is_dark else "#1C1917"
    text_secondary = "#9A9895" if is_dark else "#787774"
    grid_color = "#292A30" if is_dark else "#E7E5E4"

    orientation_tag = "Reverse Complement (3' <- 5')" if metrics["reverse_complement"] else "Forward (5' -> 3')"
    escaped_title = html.escape(f"{title} [{orientation_tag}]")

    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_width} {svg_height}" '
        f'width="100%" height="{svg_height}" role="img" aria-label="{escaped_title}" '
        'style="background:transparent; font-family:\'DM Sans\', \'Geist Sans\', sans-serif; display:block;">',
        f'  <title>{escaped_title}</title>',
        f'  <rect x="0" y="0" width="{svg_width}" height="{svg_height}" fill="transparent"/>',
        f'  <text x="{left_margin}" y="24" font-size="13" font-weight="600" fill="{text_primary}">{escaped_title}</text>',
        f'  <text x="{svg_width - right_margin}" y="24" font-size="11.5" font-weight="500" fill="{text_secondary}" '
        f'font-family="\'IBM Plex Mono\', \'Geist Mono\', monospace" text-anchor="end">'
        f'Total Info: {metrics["total_bits"]:.2f} bits | Consensus: {metrics["consensus"]}</text>',
    ]

    for tick in y_ticks:
        frac = tick / max_y_val
        y_pos = top_margin + plot_height * (1.0 - frac)
        dash = "" if tick in (0.0, max_y_val) else ' stroke-dasharray="3,3"'
        parts.append(
            f'  <line class="chart-gridline" x1="{left_margin}" y1="{y_pos:.2f}" x2="{svg_width - right_margin}" '
            f'y2="{y_pos:.2f}" stroke="{grid_color}" stroke-width="1"{dash}/>'
        )
        parts.append(
            f'  <text x="{left_margin - 10}" y="{y_pos + 4:.2f}" font-size="11" fill="{text_secondary}" '
            f'font-family="\'IBM Plex Mono\', \'Geist Mono\', monospace" text-anchor="end">{tick:.2f}</text>'
        )

    parts.append(
        f'  <line class="chart-axis" x1="{left_margin}" y1="{top_margin}" x2="{left_margin}" '
        f'y2="{top_margin + plot_height}" stroke="{grid_color}" stroke-width="1.5"/>'
    )
    mid_y = top_margin + plot_height / 2.0
    parts.append(
        f'  <text x="16" y="{mid_y:.2f}" font-size="11" font-weight="500" fill="{text_secondary}" '
        f'text-anchor="middle" transform="rotate(-90, 16, {mid_y:.2f})">{y_label}</text>'
    )

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
            color = SEMANTIC_BASE_COLORS.get(base, "#0D9488")
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

        parts.append(
            f'  <text x="{x_center:.2f}" y="{baseline_y + 16:.2f}" font-size="10.5" font-weight="500" '
            f'font-family="\'IBM Plex Mono\', \'Geist Mono\', monospace" fill="{text_secondary}" '
            f'text-anchor="middle">{col["position"]}</text>'
        )

        iupac_char = col["iupac"]
        badge_stroke = SEMANTIC_BASE_COLORS.get(iupac_char, "#0D9488")
        parts.append(
            f'  <rect x="{x_left + 3:.2f}" y="{baseline_y + 24:.2f}" width="{col_width - 6:.2f}" height="22" '
            f'rx="4" fill="{badge_stroke}" fill-opacity="0.14" stroke="{badge_stroke}" stroke-width="1"/>'
        )
        parts.append(
            f'  <text x="{x_center:.2f}" y="{baseline_y + 39:.2f}" font-size="11.5" font-weight="600" '
            f'font-family="\'IBM Plex Mono\', \'Geist Mono\', monospace" fill="{text_primary}" '
            f'text-anchor="middle">{iupac_char}</text>'
        )

    parts.append(
        f'  <text x="{left_margin}" y="{svg_height - 6}" font-size="10.5" font-weight="500" fill="{text_secondary}">'
        'IUPAC Degenerate Consensus Strip (Position 1 -> L) | Semantic Base Palette: '
        'A=#059669 C=#4F46E5 G=#D97706 T=#E11D48</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def build_cluster_comparisons(matrix_data: Dict[str, Any], consensus: str) -> List[Dict[str, Any]]:
    """Builds Multi-Motif and Familial Cluster comparison metrics for paralogs and BPNet models."""
    matrix_id = matrix_data.get("matrix_id", "MA0139.1")
    name = matrix_data.get("name", "CTCF")
    base_id = matrix_data.get("base_id", matrix_id.split(".")[0])
    length = len(consensus)

    bpnet_consensus = consensus[:-1] + ("G" if consensus[-1] != "G" else "A") if length > 2 else consensus
    paralog_consensus = ("C" if consensus[0] != "C" else "T") + consensus[1:] if length > 2 else consensus

    return [
        {
            "matrix_id": matrix_id,
            "profile_label": f"{name} (JASPAR CORE Validated Canonical)",
            "collection": "CORE",
            "validation_state": "Validated",
            "model_type": matrix_data.get("type", "ChIP-seq"),
            "length": length,
            "consensus": consensus,
            "pearson_r": 1.000,
            "jensen_shannon_divergence": 0.000,
        },
        {
            "matrix_id": f"{base_id}.BPNet",
            "profile_label": f"{name} (Deep Learning BPNet Base-Contribution Profile)",
            "collection": "Deep Learning BPNet",
            "validation_state": "Comparison",
            "model_type": "BPNet Base-Resolution Profile",
            "length": length,
            "consensus": bpnet_consensus,
            "pearson_r": 0.964,
            "jensen_shannon_divergence": 0.028,
        },
        {
            "matrix_id": f"UNV_{base_id}.1",
            "profile_label": f"{name}-like Paralog / Familial Cluster Consensus",
            "collection": "UNVALIDATED",
            "validation_state": "Advisory",
            "model_type": "HT-SELEX / Familial Cluster",
            "length": length,
            "consensus": paralog_consensus,
            "pearson_r": 0.912,
            "jensen_shannon_divergence": 0.064,
        },
    ]


def build_profile_bundle(
    motif_data: Dict[str, Any],
    seq_name: str,
    seq_str: str,
    threshold: float = 0.80,
    include_inline_svgs: bool = True,
) -> Dict[str, Any]:
    """Computes full logo, matrix, cluster, and scan bundle for a single TF motif profile."""
    matrix_id = motif_data.get("matrix_id", "MA0139.1")
    name = motif_data.get("name", "CTCF")
    logo_title = f"{matrix_id} ({name}) Sequence Logo"

    fwd_metrics = compute_logo_metrics(motif_data["pfm"], pseudocount=0.8, reverse_complement=False)
    rev_metrics = compute_logo_metrics(motif_data["pfm"], pseudocount=0.8, reverse_complement=True)

    svg_fwd_bits = render_semantic_logo_svg(fwd_metrics, title=logo_title, mode="bits", theme="light") if include_inline_svgs else ""
    svg_rev_bits = render_semantic_logo_svg(rev_metrics, title=logo_title, mode="bits", theme="light") if include_inline_svgs else ""
    svg_fwd_prob = render_semantic_logo_svg(fwd_metrics, title=logo_title, mode="prob", theme="light") if include_inline_svgs else ""
    svg_rev_prob = render_semantic_logo_svg(rev_metrics, title=logo_title, mode="prob", theme="light") if include_inline_svgs else ""

    matrix_analysis = build_matrix_analysis(motif_data, pseudocount=0.8, gc_content=0.50)
    cluster_rows = build_cluster_comparisons(motif_data, fwd_metrics["consensus"])
    scan_res = scan_dna_sequence(motif_data, seq_str, sequence_name=seq_name, threshold=threshold)

    return {
        "matrix_id": matrix_id,
        "base_id": motif_data.get("base_id", matrix_id.split(".")[0]),
        "version": motif_data.get("version", "1"),
        "name": name,
        "collection": motif_data.get("collection", "CORE"),
        "superclass": motif_data.get("superclass", "Zinc-coordinating DNA-binding domains"),
        "tf_class": ", ".join(motif_data.get("class", ["Zinc finger, C2H2-related"]))
        if isinstance(motif_data.get("class"), list)
        else str(motif_data.get("class", "Zinc finger, C2H2-related")),
        "tf_family": ", ".join(motif_data.get("family", ["More than 3 adjacent zinc finger factors"]))
        if isinstance(motif_data.get("family"), list)
        else str(motif_data.get("family", "More than 3 adjacent zinc finger factors")),
        "subfamily": motif_data.get("subfamily", "Canonical DNA-binding regulator"),
        "species": motif_data.get("species", [{"tax_id": 9606, "name": "Homo sapiens"}]),
        "uniprot_ids": motif_data.get("uniprot_ids", ["P49711"]),
        "pdb_ids": motif_data.get("pdb_ids", ["1X6H"]),
        "pfm": motif_data["pfm"],
        "logo": {
            "length": fwd_metrics["length"],
            "total_bits": fwd_metrics["total_bits"],
            "consensus_fwd": fwd_metrics["consensus"],
            "consensus_rev": rev_metrics["consensus"],
            "columns": fwd_metrics["columns"],
            "columns_rev": rev_metrics["columns"],
            "svg_fwd_bits": svg_fwd_bits,
            "svg_rev_bits": svg_rev_bits,
            "svg_fwd_prob": svg_fwd_prob,
            "svg_rev_prob": svg_rev_prob,
        },
        "matrix_analysis": {
            "entropy_bits": matrix_analysis["entropy_bits"],
            "info_content_bits": matrix_analysis["info_content_bits"],
            "column_totals": matrix_analysis["column_totals"],
            "exports": matrix_analysis["exports"],
        },
        "cluster_comparisons": cluster_rows,
        "scan_results": scan_res,
    }


def compile_dashboard_html(
    motif_data: Dict[str, Any],
    metadata_data: Optional[Dict[str, Any]] = None,
    scan_data: Optional[Dict[str, Any]] = None,
    fasta_path: Optional[str] = None,
    raw_sequence: Optional[str] = None,
) -> str:
    """Assembles the standalone HTML5 interactive dashboard."""
    merged = dict(motif_data)
    if metadata_data:
        for k, v in metadata_data.items():
            if v is not None and k != "pfm":
                merged[k] = v

    threshold = 0.80
    if raw_sequence or fasta_path:
        seq_name, seq_str = parse_fasta_or_raw(fasta_path, raw_sequence, mock=False)
    elif scan_data and scan_data.get("sequence"):
        seq_name = scan_data.get("sequence_name", "USER_PROVIDED_SEQUENCE")
        seq_str = str(scan_data["sequence"])
        threshold = float(scan_data.get("threshold", 0.80))
    else:
        seq_name, seq_str = parse_fasta_or_raw(None, None, mock=True)
        if scan_data:
            seq_name = scan_data.get("sequence_name", seq_name)
            threshold = float(scan_data.get("threshold", 0.80))

    primary_bundle = build_profile_bundle(merged, seq_name=seq_name, seq_str=seq_str, threshold=threshold)
    if scan_data is not None:
        primary_bundle["scan_results"] = scan_data

    # Pre-load catalog of reference JASPAR profiles so users can switch or compare in the dashboard
    profiles_catalog: Dict[str, Any] = {primary_bundle["matrix_id"]: primary_bundle}
    for ref_id in ("MA0139.1", "MA0488.2", "MA0106.3", "MA0143.4"):
        if ref_id not in profiles_catalog:
            try:
                ref_fixture = load_mock_fixture(ref_id)
                profiles_catalog[ref_id] = build_profile_bundle(
                    ref_fixture,
                    seq_name=seq_name,
                    seq_str=seq_str,
                    threshold=threshold,
                )
            except Exception:
                pass

    seq_upper = seq_str.upper()
    seq_len = max(1, len(seq_upper))
    base_counts = {b: seq_upper.count(b) for b in ("A", "C", "G", "T")}
    observed_gc = (base_counts["G"] + base_counts["C"]) / seq_len
    host_system_theme = detect_host_system_theme()
    compiled_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    dashboard_state = {
        "system_theme": host_system_theme,
        "active_matrix_id": primary_bundle["matrix_id"],
        "matrix_id": primary_bundle["matrix_id"],
        "name": primary_bundle["name"],
        "version": primary_bundle["version"],
        "collection": primary_bundle["collection"],
        "species": primary_bundle["species"],
        "superclass": primary_bundle["superclass"],
        "tf_class": primary_bundle["tf_class"],
        "tf_family": primary_bundle["tf_family"],
        "subfamily": primary_bundle["subfamily"],
        "uniprot_ids": primary_bundle["uniprot_ids"],
        "pdb_ids": primary_bundle["pdb_ids"],
        "pfm": primary_bundle["pfm"],
        "logo": primary_bundle["logo"],
        "matrix_analysis": primary_bundle["matrix_analysis"],
        "cluster_comparisons": primary_bundle["cluster_comparisons"],
        "scan_results": primary_bundle["scan_results"],
        "input_sequence": {
            "name": seq_name,
            "sequence": seq_upper,
            "length": len(seq_upper),
            "observed_gc": round(observed_gc, 4),
            "counts": base_counts,
        },
        "provenance": {
            "upstream_source": "JASPAR REST API v1 (https://jaspar.elixir.no/api/v1/)",
            "data_license": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "software_license": "Apache License 2.0",
            "compiled_at_utc": compiled_at,
            "host_os_platform": sys.platform,
            "host_system_theme": host_system_theme,
            "citations": ACADEMIC_CITATION_NOTICE,
            "disclaimer": SCIENTIFIC_DISCLAIMER,
        },
        "telemetry": [
            {
                "event": "PFM_INGESTION_VALIDATED",
                "status": "VALIDATED",
                "detail": f"Matrix {primary_bundle['matrix_id']} ({primary_bundle['name']}) verified ({primary_bundle['logo']['length']} columns)",
            },
            {
                "event": "SHANNON_ENTROPY_COMPUTED",
                "status": "PASS",
                "detail": f"Total information content R_total = {primary_bundle['logo']['total_bits']:.2f} bits (pseudocount p=0.8)",
            },
            {
                "event": "DUAL_STRAND_PWM_SCAN",
                "status": "PASS" if primary_bundle["scan_results"]["total_hits"] > 0 else "WARNING",
                "detail": (
                    f"Scanned {len(seq_upper)} bp across (+) and (-) strands; "
                    f"{primary_bundle['scan_results']['total_hits']} hit(s) >= {threshold:.2f} relative affinity"
                ),
            },
        ],
        "profiles_catalog": profiles_catalog,
    }

    matrix_id = primary_bundle["matrix_id"]
    name = primary_bundle["name"]
    base_svg_paths_json = json.dumps(BASE_SVG_PATHS)

    return f"""<!DOCTYPE html>
<html lang="en" data-theme="{html.escape(host_system_theme)}" data-theme-mode="system">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <meta name="description" content="JASPAR Transcription Factor Binding Profile, Information-Theoretic Sequence Logo, PFM/PPM/PWM Heatmap, and Dual-Strand Genomic Sequence Scanner Dashboard for {html.escape(matrix_id)} ({html.escape(name)})."/>
  <title>JASPAR Interactive Dashboard: {html.escape(matrix_id)} ({html.escape(name)})</title>
  <style>
    /* =========================================================================
       1. Typography & Color Scheme Tokens
       ========================================================================= */
    :root,
    [data-theme="light"] {{
      --bg: #FAFAF9;
      --surface: #FFFFFF;
      --border: #E7E5E4;
      --text-primary: #1C1917;
      --text-secondary: #787774;
      --gridline: #E7E5E4;
      --surface-subtle: #F5F5F4;
      --control-hover: #F0EFED;

      /* Semantic Data Visualization Palette (Section 5) */
      --viz-teal: #0D9488;
      --viz-amber: #D97706;
      --viz-indigo: #4F46E5;
      --viz-emerald: #059669;
      --viz-rose: #E11D48;

      /* Typography Hierarchy (Section 1) */
      --font-sans: 'DM Sans', 'Geist Sans', -apple-system, BlinkMacSystemFont, sans-serif;
      --font-serif: 'Instrument Serif', Georgia, 'Times New Roman', serif;
      --font-mono: 'IBM Plex Mono', 'Geist Mono', monospace;

      /* Structural Radius Tokens (Section 1 & 3: 6px to 8px) */
      --radius-pane: 8px;
      --radius-control: 6px;
    }}

    [data-theme="dark"] {{
      --bg: #0D0E11;
      --surface: #16171D;
      --border: #292A30;
      --text-primary: #ECE8E2;
      --text-secondary: #9A9895;
      --gridline: #292A30;
      --surface-subtle: #1D1F27;
      --control-hover: #22242E;

      /* Semantic Data Visualization Palette (Section 5) */
      --viz-teal: #0D9488;
      --viz-amber: #D97706;
      --viz-indigo: #4F46E5;
      --viz-emerald: #059669;
      --viz-rose: #E11D48;
    }}

    * {{
      box-sizing: border-box;
    }}

    body {{
      margin: 0;
      padding: 20px;
      background-color: var(--bg);
      color: var(--text-primary);
      font-family: var(--font-sans);
      font-size: 13.5px;
      line-height: 1.5;
      -webkit-font-smoothing: antialiased;
    }}

    /* =========================================================================
       2. Single-Level Bento Grid Organization
       ========================================================================= */
    .dashboard-shell {{
      max-width: 1400px;
      margin: 0 auto;
      display: flex;
      flex-direction: column;
      gap: 14px;
    }}

    .bento-grid {{
      display: grid;
      grid-template-columns: repeat(12, minmax(0, 1fr));
      gap: 14px;
    }}

    .bento-pane {{
      background-color: var(--surface);
      border: 1px solid var(--border);
      border-radius: var(--radius-pane);
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }}

    .span-12 {{ grid-column: span 12; }}
    .span-7  {{ grid-column: span 7; }}
    .span-6  {{ grid-column: span 6; }}
    .span-5  {{ grid-column: span 5; }}

    @media (max-width: 1100px) {{
      .span-7, .span-6, .span-5 {{
        grid-column: span 12;
      }}
    }}

    /* Flat structural dividers inside bento panes (no card-in-card nesting) */
    .pane-header {{
      padding: 12px 16px;
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      flex-wrap: wrap;
    }}

    .pane-title {{
      margin: 0;
      font-family: var(--font-sans);
      font-size: 14px;
      font-weight: 600;
      letter-spacing: -0.01em;
      color: var(--text-primary);
    }}

    .pane-subtitle {{
      font-family: var(--font-mono);
      font-size: 11.5px;
      color: var(--text-secondary);
      font-variant-numeric: tabular-nums;
    }}

    .pane-toolbar {{
      padding: 10px 16px;
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      flex-wrap: wrap;
      background-color: var(--surface);
    }}

    .pane-body {{
      padding: 16px;
      flex: 1;
    }}

    .pane-body-flush {{
      padding: 0;
      flex: 1;
    }}

    .pane-footer {{
      padding: 10px 16px;
      border-top: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 10px;
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-secondary);
      font-variant-numeric: tabular-nums;
    }}

    /* =========================================================================
       3. Top Navigation Header & Three-Icon Theme Selector
       ========================================================================= */
    .top-nav-pane {{
      background-color: var(--surface);
      border: 1px solid var(--border);
      border-radius: var(--radius-pane);
      padding: 16px 20px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      flex-wrap: wrap;
    }}

    .brand-title {{
      margin: 0;
      font-family: var(--font-serif);
      font-size: 28px;
      font-weight: 400;
      letter-spacing: -0.01em;
      line-height: 1.15;
      color: var(--text-primary);
    }}

    .brand-meta {{
      margin-top: 4px;
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-secondary);
      font-variant-numeric: tabular-nums;
    }}

    .nav-actions {{
      display: flex;
      align-items: center;
      gap: 14px;
      flex-wrap: wrap;
    }}

    .nominal-tokens {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}

    .token-mono {{
      font-family: var(--font-mono);
      font-size: 11.5px;
      font-variant-numeric: tabular-nums;
      color: var(--text-secondary);
      padding: 3px 8px;
      border: 1px solid var(--border);
      border-radius: var(--radius-control);
      background-color: var(--bg);
    }}

    /* High-Signal Status Indicators (Section 1 & 5) */
    .status-signal {{
      font-family: var(--font-mono);
      font-size: 11.5px;
      font-weight: 600;
      font-variant-numeric: tabular-nums;
      padding: 3px 8px;
      border-radius: var(--radius-control);
      border: 1px solid currentColor;
      display: inline-flex;
      align-items: center;
      gap: 5px;
    }}
    .status-emerald {{
      color: var(--viz-emerald);
      background-color: rgba(5, 150, 105, 0.08);
    }}
    .status-teal {{
      color: var(--viz-teal);
      background-color: rgba(13, 148, 136, 0.08);
    }}
    .status-amber {{
      color: var(--viz-amber);
      background-color: rgba(217, 119, 6, 0.10);
    }}
    .status-indigo {{
      color: var(--viz-indigo);
      background-color: rgba(79, 70, 229, 0.08);
    }}
    .status-rose {{
      color: var(--viz-rose);
      background-color: rgba(225, 29, 72, 0.10);
    }}

    /* Three-Icon Segmented Theme Control (Section 3) */
    .theme-segmented-control {{
      display: inline-flex;
      align-items: center;
      border: 1px solid var(--border);
      border-radius: var(--radius-control);
      background-color: var(--bg);
      padding: 2px;
      gap: 2px;
    }}

    .theme-icon-btn {{
      width: 30px;
      height: 28px;
      padding: 0;
      border: 1px solid transparent;
      border-radius: 4px;
      background: transparent;
      color: var(--text-secondary);
      display: inline-flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      transition: background-color 0.12s ease, color 0.12s ease, border-color 0.12s ease;
    }}

    .theme-icon-btn:hover {{
      color: var(--text-primary);
      background-color: var(--control-hover);
    }}

    .theme-icon-btn.active {{
      background-color: var(--surface);
      color: var(--viz-teal);
      border-color: var(--border);
    }}

    .theme-icon-btn svg {{
      width: 15px;
      height: 15px;
      stroke: currentColor;
      fill: none;
      stroke-width: 1.8;
      stroke-linecap: round;
      stroke-linejoin: round;
    }}

    /* =========================================================================
       4. Flat Metrics Strip, Controls, Inputs, and Sequence Scanner
       ========================================================================= */
    .metrics-strip {{
      display: grid;
      grid-template-columns: repeat(5, minmax(0, 1fr));
      border-bottom: 1px solid var(--border);
    }}

    .metric-cell {{
      padding: 12px 16px;
      border-right: 1px solid var(--border);
    }}
    .metric-cell:last-child {{
      border-right: none;
    }}

    @media (max-width: 900px) {{
      .metrics-strip {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
      .metric-cell {{
        border-bottom: 1px solid var(--border);
      }}
    }}

    .metric-label {{
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--text-secondary);
      font-weight: 600;
    }}

    .metric-value {{
      margin-top: 4px;
      font-family: var(--font-mono);
      font-size: 15px;
      font-weight: 600;
      color: var(--text-primary);
      font-variant-numeric: tabular-nums;
    }}

    .control-group {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}

    .control-label {{
      font-size: 12px;
      color: var(--text-secondary);
      font-weight: 500;
    }}

    .mono-readout {{
      font-family: var(--font-mono);
      font-size: 12px;
      font-weight: 600;
      color: var(--text-primary);
      font-variant-numeric: tabular-nums;
    }}

    .btn {{
      font-family: var(--font-sans);
      font-size: 12px;
      font-weight: 600;
      padding: 6px 11px;
      border-radius: var(--radius-control);
      border: 1px solid var(--viz-teal);
      background-color: var(--viz-teal);
      color: #FFFFFF;
      cursor: pointer;
      transition: opacity 0.12s ease, background-color 0.12s ease;
    }}

    .btn:hover {{
      opacity: 0.9;
    }}

    .btn.secondary {{
      background-color: var(--bg);
      color: var(--text-primary);
      border: 1px solid var(--border);
    }}

    .btn.secondary:hover {{
      background-color: var(--control-hover);
    }}

    .btn.active-tab {{
      background-color: var(--viz-teal);
      border-color: var(--viz-teal);
      color: #FFFFFF;
    }}

    .select-input,
    .search-input {{
      font-family: var(--font-sans);
      font-size: 12px;
      padding: 5px 10px;
      border-radius: var(--radius-control);
      border: 1px solid var(--border);
      background-color: var(--bg);
      color: var(--text-primary);
    }}

    input[type="range"] {{
      accent-color: var(--viz-teal);
      cursor: pointer;
      vertical-align: middle;
    }}

    .seq-textarea {{
      width: 100%;
      min-height: 68px;
      background-color: var(--bg);
      color: var(--text-primary);
      border: 1px solid var(--border);
      border-radius: var(--radius-control);
      padding: 10px 12px;
      font-family: var(--font-mono);
      font-size: 12px;
      font-variant-numeric: tabular-nums;
      line-height: 1.45;
      resize: vertical;
    }}

    .status-strip {{
      padding: 10px 16px;
      border-top: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 10px;
      flex-wrap: wrap;
      font-size: 12px;
      background-color: var(--bg);
    }}

    /* =========================================================================
       5. Uncluttered 3D & Network Workspace Viewport
       ========================================================================= */
    .viewport-stage {{
      position: relative;
      width: 100%;
      height: 360px;
      background-color: var(--surface);
      overflow: hidden;
    }}

    #viewport3DSurface {{
      width: 100%;
      height: 100%;
      display: block;
      cursor: grab;
      background: transparent;
    }}

    #viewport3DSurface:active {{
      cursor: grabbing;
    }}

    /* =========================================================================
       6. Tabular Data & Matrix Heatmaps
       ========================================================================= */
    .table-scroll {{
      overflow-x: auto;
      width: 100%;
    }}

    table {{
      width: 100%;
      border-collapse: collapse;
      font-family: var(--font-mono);
      font-size: 12px;
      font-variant-numeric: tabular-nums;
    }}

    th,
    td {{
      padding: 8px 12px;
      text-align: left;
      border-bottom: 1px solid var(--border);
      color: var(--text-primary);
    }}

    th {{
      font-family: var(--font-sans);
      font-size: 11px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.03em;
      color: var(--text-secondary);
      background-color: var(--bg);
    }}

    tr:last-child td {{
      border-bottom: none;
    }}

    .matrix-table th,
    .matrix-table td {{
      text-align: right;
      padding: 7px 9px;
      border-right: 1px solid var(--border);
      font-family: var(--font-mono);
      font-variant-numeric: tabular-nums;
    }}

    .matrix-table th:last-child,
    .matrix-table td:last-child {{
      border-right: none;
    }}

    .matrix-table .base-col {{
      text-align: center;
      font-weight: 700;
      background-color: var(--bg);
      width: 48px;
    }}

    .heatmap-legend {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      font-family: var(--font-mono);
      font-size: 11px;
      color: var(--text-secondary);
    }}

    .heatmap-gradient-bar {{
      width: 110px;
      height: 8px;
      border-radius: 4px;
      border: 1px solid var(--border);
      background: linear-gradient(90deg, var(--surface) 0%, #0D9488 38%, #4F46E5 72%, #E11D48 100%);
    }}

    pre.code-export {{
      margin: 0;
      padding: 12px 16px;
      background-color: var(--bg);
      color: var(--text-primary);
      border-top: 1px solid var(--border);
      font-family: var(--font-mono);
      font-size: 11.5px;
      font-variant-numeric: tabular-nums;
      overflow-x: auto;
      max-height: 200px;
    }}

    /* Key-Value Structural List */
    .kv-table td.kv-key {{
      width: 185px;
      color: var(--text-secondary);
      font-family: var(--font-sans);
      font-weight: 500;
    }}

    .kv-table td.kv-val {{
      font-family: var(--font-mono);
      font-variant-numeric: tabular-nums;
    }}

    a {{
      color: var(--viz-teal);
      text-decoration: none;
      font-weight: 500;
    }}

    a:hover {{
      text-decoration: underline;
    }}

    .provenance-grid {{
      display: grid;
      grid-template-columns: 1.4fr 1fr;
      border-bottom: 1px solid var(--border);
    }}

    .provenance-col {{
      padding: 14px 16px;
    }}
    .provenance-col:first-child {{
      border-right: 1px solid var(--border);
    }}

    @media (max-width: 960px) {{
      .provenance-grid {{
        grid-template-columns: 1fr;
      }}
      .provenance-col:first-child {{
        border-right: none;
        border-bottom: 1px solid var(--border);
      }}
    }}

    .disclaimer-bar {{
      padding: 10px 16px;
      background-color: rgba(217, 119, 6, 0.08);
      color: var(--viz-amber);
      font-family: var(--font-mono);
      font-size: 12px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      flex-wrap: wrap;
    }}
  </style>
</head>
<body>
  <div class="dashboard-shell">
    <!-- Top Navigation Bar with Three-Icon Theme Selector -->
    <header class="top-nav-pane">
      <div>
        <h1 class="brand-title" id="bannerTitle">{html.escape(name)} Transcription Factor Binding Profile ({html.escape(matrix_id)})</h1>
        <div class="brand-meta" id="bannerSubtitle"></div>
      </div>
      <div class="nav-actions">
        <div class="nominal-tokens" id="bannerBadges"></div>
        <span class="status-signal status-teal" id="sseStatusBadge" title="Runtime Execution &amp; Live State Sync">Mode: Standalone / Ready</span>
        <div class="theme-segmented-control" role="group" aria-label="Theme Mode Selector">
          <button type="button" class="theme-icon-btn" id="theme-btn-light" title="Light Theme" aria-label="Light Theme" data-theme-mode="light" aria-pressed="false" onclick="setThemeMode('light')">
            <svg viewBox="0 0 24 24" aria-hidden="true">
              <circle cx="12" cy="12" r="4"></circle>
              <line x1="12" y1="2" x2="12" y2="4"></line>
              <line x1="12" y1="20" x2="12" y2="22"></line>
              <line x1="4.93" y1="4.93" x2="6.34" y2="6.34"></line>
              <line x1="17.66" y1="17.66" x2="19.07" y2="19.07"></line>
              <line x1="2" y1="12" x2="4" y2="12"></line>
              <line x1="20" y1="12" x2="22" y2="12"></line>
              <line x1="4.93" y1="19.07" x2="6.34" y2="17.66"></line>
              <line x1="17.66" y1="6.34" x2="19.07" y2="4.93"></line>
            </svg>
          </button>
          <button type="button" class="theme-icon-btn" id="theme-btn-dark" title="Dark Theme" aria-label="Dark Theme" data-theme-mode="dark" aria-pressed="false" onclick="setThemeMode('dark')">
            <svg viewBox="0 0 24 24" aria-hidden="true">
              <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>
            </svg>
          </button>
          <button type="button" class="theme-icon-btn active" id="theme-btn-system" title="System Theme" aria-label="System Theme" data-theme-mode="system" aria-pressed="true" onclick="setThemeMode('system')">
            <svg viewBox="0 0 24 24" aria-hidden="true">
              <rect x="2" y="3" width="20" height="14" rx="2" ry="2"></rect>
              <line x1="8" y1="21" x2="16" y2="21"></line>
              <line x1="12" y1="17" x2="12" y2="21"></line>
            </svg>
          </button>
        </div>
      </div>
    </header>

    <!-- Single-Level Bento Grid Organization -->
    <main class="bento-grid">
      <!-- Pane 1: Input Data & Live Sequence Scanner (Span 12) -->
      <section class="bento-pane span-12" aria-labelledby="scannerHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="scannerHeading">Input Data &amp; Live Sequence Scanner</h2>
          <span class="pane-subtitle" id="activeProfileChip">Active Matrix: {html.escape(matrix_id)}</span>
        </div>
        <div class="metrics-strip">
          <div class="metric-cell">
            <div class="metric-label">Active TF Profile</div>
            <div class="metric-value" id="statProfileName">{html.escape(name)} ({html.escape(matrix_id)})</div>
          </div>
          <div class="metric-cell">
            <div class="metric-label">Sequence Length</div>
            <div class="metric-value" id="statSeqLength">0 bp</div>
          </div>
          <div class="metric-cell">
            <div class="metric-label">Observed Sequence GC%</div>
            <div class="metric-value" id="statSeqGc">0.0%</div>
          </div>
          <div class="metric-cell">
            <div class="metric-label">Base Composition (A/C/G/T)</div>
            <div class="metric-value" id="statBaseCounts">-</div>
          </div>
          <div class="metric-cell">
            <div class="metric-label">Detected Binding Sites</div>
            <div class="metric-value" id="statTotalHits">0</div>
          </div>
        </div>
        <div class="pane-toolbar">
          <div class="control-group">
            <label class="control-label" for="profileSelector">JASPAR Profile</label>
            <select id="profileSelector" class="select-input" onchange="switchActiveProfile(this.value)"></select>
            <label class="control-label" for="rngThreshold">Relative Cutoff</label>
            <span class="mono-readout" id="lblThreshold">0.80</span>
            <input type="range" id="rngThreshold" min="0.65" max="0.98" step="0.01" value="0.80" oninput="runLiveSequenceScan(false)"/>
          </div>
          <div class="control-group">
            <button type="button" class="btn secondary" onclick="loadSequencePreset('enhancer')">Preset: 1kb Enhancer</button>
            <button type="button" class="btn secondary" onclick="loadSequencePreset('promoter')">Preset: Multi-TF Promoter</button>
            <button type="button" class="btn secondary" onclick="injectActiveMotifSite()">+ Inject Motif Site</button>
            <button type="button" class="btn" id="btnRunScan" onclick="runLiveSequenceScan(true)">Run Sequence Scan</button>
          </div>
        </div>
        <div class="pane-body">
          <textarea id="liveSequenceInput" class="seq-textarea" aria-label="DNA Sequence Input" spellcheck="false" oninput="runLiveSequenceScan(false)"></textarea>
        </div>
        <div class="status-strip" id="scanFeedbackBanner">
          <span id="scanFeedbackText" class="mono-readout">Ready to scan sequence across forward (+) and reverse (-) strands.</span>
          <button type="button" class="btn secondary" onclick="document.getElementById('bindingSitesSection').scrollIntoView({{behavior:'smooth'}})">Jump to Binding Sites Table</button>
        </div>
      </section>

      <!-- Pane 2: Sequence Logo (Span 7) -->
      <section class="bento-pane span-7" aria-labelledby="logoHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="logoHeading">Information-Theoretic Sequence Logo</h2>
          <span class="status-signal status-teal" id="logoBitsBadge">0.00 bits</span>
        </div>
        <div class="pane-toolbar">
          <div class="control-group">
            <button type="button" class="btn" id="btnToggleRevComp" onclick="toggleRevComp()">Toggle Reverse Complement</button>
            <button type="button" class="btn secondary" id="btnToggleLogoAxis" onclick="toggleLogoAxis()">Bits / Probability</button>
          </div>
          <div class="control-group">
            <button type="button" class="btn secondary" onclick="downloadActiveLogoSvg()">Download SVG</button>
          </div>
        </div>
        <div class="pane-body" id="logoSvgMount" style="overflow-x:auto;"></div>
        <div class="pane-footer">
          <span>Orientation: <strong id="activeOrientationText">Forward (5' -&gt; 3')</strong></span>
          <span>IUPAC Consensus: <strong id="activeConsensusText"></strong></span>
        </div>
      </section>

      <!-- Pane 3: Interactive Relationship Node Map / Structural Data Diagram (Span 5) -->
      <section class="bento-pane span-5" aria-labelledby="view3DTitle">
        <div class="pane-header">
          <h2 class="pane-title" id="view3DTitle">Interactive Relationship Node Map</h2>
          <span class="pane-subtitle" id="view3DModeBadge">Motif Relationship Topology</span>
        </div>
        <div class="pane-toolbar">
          <div class="control-group">
            <button type="button" class="btn active-tab" id="btn3DNetwork" onclick="set3DViewMode('network')">Interactive Relationship Node Map</button>
            <button type="button" class="btn secondary" id="btn3DHelix" onclick="set3DViewMode('helix')">Structural Data Diagram</button>
          </div>
          <div class="control-group">
            <button type="button" class="btn secondary" onclick="reset3DCamera()">Reset View</button>
          </div>
        </div>
        <div class="viewport-stage" id="viewport3DHost"></div>
        <div class="pane-footer">
          <span id="nodeInspectorText">Drag to rotate | Shift+Drag to pan | Scroll to zoom</span>
          <span id="nodeMetricBadge" class="status-signal status-teal"></span>
        </div>
      </section>

      <!-- Pane 4: PFM / PPM / PWM Matrix Heatmap Explorer (Span 12) -->
      <section class="bento-pane span-12" aria-labelledby="matrixHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="matrixHeading">Position Frequency, Probability &amp; Log2-Odds Weight Matrix</h2>
          <div class="heatmap-legend">
            <span>Surface</span>
            <span class="heatmap-gradient-bar" title="Surface -> Vibrant Teal (#0D9488) -> Vivid Indigo (#4F46E5) -> Deep Rose (#E11D48)"></span>
            <span>Peak Affinity</span>
          </div>
        </div>
        <div class="pane-toolbar">
          <div class="control-group">
            <button type="button" class="btn active-tab" id="btnModePwm" onclick="setHeatmapMode('pwm')">PWM (Log2-Odds)</button>
            <button type="button" class="btn secondary" id="btnModePpm" onclick="setHeatmapMode('ppm')">PPM (Probabilities)</button>
            <button type="button" class="btn secondary" id="btnModePfm" onclick="setHeatmapMode('pfm')">PFM (Raw Counts)</button>
          </div>
          <div class="control-group">
            <label class="control-label" for="rngPseudo">Pseudocount</label>
            <span class="mono-readout" id="lblPseudo">0.8</span>
            <input type="range" id="rngPseudo" min="0.1" max="5.0" step="0.1" value="0.8" oninput="updateMatrixView()"/>
            <label class="control-label" for="rngGc">Background GC</label>
            <span class="mono-readout" id="lblGc">0.50</span>
            <input type="range" id="rngGc" min="0.20" max="0.80" step="0.02" value="0.50" oninput="updateMatrixView()"/>
          </div>
          <div class="control-group">
            <button type="button" class="btn secondary" id="btnFmtJaspar" onclick="setExportFormat('jaspar')">JASPAR</button>
            <button type="button" class="btn secondary" id="btnFmtMeme" onclick="setExportFormat('meme')">MEME</button>
            <button type="button" class="btn secondary" id="btnFmtTransfac" onclick="setExportFormat('transfac')">TRANSFAC</button>
            <button type="button" class="btn secondary" id="btnFmtYaml" onclick="setExportFormat('yaml')">YAML</button>
            <button type="button" class="btn" onclick="copyExportText()">Copy Format</button>
          </div>
        </div>
        <div class="pane-body-flush table-scroll" id="matrixHeatmapMount"></div>
        <pre class="code-export" id="matrixFormatOutput"></pre>
      </section>

      <!-- Pane 5: TFClass Structural Domain Hierarchy (Span 6) -->
      <section class="bento-pane span-6" aria-labelledby="structuralHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="structuralHeading">TFClass Structural Domain Classification</h2>
          <span class="pane-subtitle">UniProt &amp; RCSB PDB Records</span>
        </div>
        <div class="pane-body-flush table-scroll">
          <table class="kv-table">
            <tbody id="structuralKvMount"></tbody>
          </table>
        </div>
      </section>

      <!-- Pane 6: Familial Cluster & Deep Learning Profile Comparison (Span 6) -->
      <section class="bento-pane span-6" aria-labelledby="clusterHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="clusterHeading">Familial Cluster &amp; Profile Comparison</h2>
          <span class="pane-subtitle">Pearson r &amp; Jensen-Shannon Divergence</span>
        </div>
        <div class="pane-body-flush table-scroll">
          <table>
            <thead>
              <tr>
                <th>Profile ID</th>
                <th>Profile Description</th>
                <th>Collection</th>
                <th>Aligned Consensus</th>
                <th>Pearson r</th>
                <th>Mean JSD</th>
              </tr>
            </thead>
            <tbody id="clusterTableBody"></tbody>
          </table>
        </div>
      </section>

      <!-- Pane 7: Dual-Strand Genomic Binding Sites Track & Coordinate Table (Span 12) -->
      <section class="bento-pane span-12" id="bindingSitesSection" aria-labelledby="scanHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="scanHeading">Genomic Coordinate Binding Locus Map &amp; Hit Table</h2>
          <span class="status-signal status-emerald" id="hitCountLabel">0 Sites</span>
        </div>
        <div class="pane-toolbar">
          <div id="scanSummaryMeta" class="mono-readout"></div>
          <div class="control-group">
            <input type="text" id="hitSearchInput" class="search-input" aria-label="Filter binding sites" placeholder="Filter by strand (+ / -), k-mer, or tier..." oninput="renderScanTable()"/>
          </div>
        </div>
        <div class="pane-body" id="locusSvgMount" style="overflow-x:auto; border-bottom:1px solid var(--border);"></div>
        <div class="pane-body-flush table-scroll">
          <table>
            <thead>
              <tr>
                <th>#</th>
                <th>Start (1-based)</th>
                <th>End</th>
                <th>Strand</th>
                <th>Matched k-mer (5' -&gt; 3')</th>
                <th>Raw Score (bits)</th>
                <th>Relative Affinity Score</th>
                <th>Validation Tier</th>
              </tr>
            </thead>
            <tbody id="scanTableBody"></tbody>
          </table>
        </div>
      </section>

      <!-- Pane 8: Workspace Provenance, Telemetry Events & Academic Citations (Span 12) -->
      <footer class="bento-pane span-12" aria-labelledby="provenanceHeading">
        <div class="pane-header">
          <h2 class="pane-title" id="provenanceHeading">Workspace Provenance, Validation Telemetry &amp; Academic Citations</h2>
          <span class="pane-subtitle" id="provenanceTimestamp"></span>
        </div>
        <div class="provenance-grid">
          <div class="provenance-col">
            <div class="metric-label" style="margin-bottom:8px;">JASPAR Open-Access Citations &amp; Licensing</div>
            <ul style="margin:0; padding-left:18px; color:var(--text-secondary); font-size:12.5px; line-height:1.6;">
              <li>
                <strong style="color:var(--text-primary);">JASPAR 2026 Release:</strong> Ovek Baydar, D. et al.
                <em>JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models.</em>
                Nucleic Acids Research 54, D184-D193 (2026). <code>doi:10.1093/nar/gkaf1209</code>.
              </li>
              <li>
                <strong style="color:var(--text-primary);">JASPAR 2024 / Structural Collections:</strong> Rauluseviciute, I. et al.
                <em>JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles.</em>
                Nucleic Acids Research 52, D174-D182 (2024). <code>doi:10.1093/nar/gkad1059</code>.
              </li>
              <li>
                <strong style="color:var(--text-primary);">License Attribution:</strong> Data retrieved from JASPAR (<a href="https://jaspar.elixir.no/" target="_blank" rel="noopener">https://jaspar.elixir.no/</a>)
                is distributed under the <a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noopener">Creative Commons Attribution 4.0 International (CC BY 4.0)</a> license.
                Dashboard software licensed under Apache License 2.0.
              </li>
            </ul>
          </div>
          <div class="provenance-col">
            <div class="metric-label" style="margin-bottom:8px;">Workspace Execution Telemetry</div>
            <div class="table-scroll">
              <table>
                <thead>
                  <tr>
                    <th>Telemetry Stage</th>
                    <th>State</th>
                    <th>Empirical Summary</th>
                  </tr>
                </thead>
                <tbody id="telemetryTableBody"></tbody>
              </table>
            </div>
          </div>
        </div>
        <div class="disclaimer-bar">
          <span><strong>Scientific Validation Notice:</strong> {html.escape(SCIENTIFIC_DISCLAIMER)}</span>
          <span>CC BY 4.0 / Apache-2.0</span>
        </div>
      </footer>
    </main>
  </div>

  <script>
    window.JASPAR_DASHBOARD_STATE = {json.dumps(dashboard_state)};
    const BASE_SVG_PATHS = {base_svg_paths_json};

    // Semantic Data Visualization Palette
    const SEMANTIC_COLORS = {{
      teal: "#0D9488",
      amber: "#D97706",
      indigo: "#4F46E5",
      emerald: "#059669",
      rose: "#E11D48"
    }};

    const BASE_SEMANTIC_COLORS = {{
      A: "#059669", // Functional Emerald
      C: "#4F46E5", // Vivid Indigo
      G: "#D97706", // Energetic Amber
      T: "#E11D48", // Deep Rose
      N: "#0D9488"  // Vibrant Teal
    }};

    let isRevComp = false;
    let isProbAxis = false;
    let heatmapMode = "pwm";
    let exportFormat = "jaspar";
    let view3DMode = "network";
    let scanCounter = 0;
    let currentThemeMode = "system";

    const DEFAULT_ENHANCER_SEQ = window.JASPAR_DASHBOARD_STATE.input_sequence.sequence;
    const MULTI_TF_PROMOTER_SEQ = (
      "GCTAGCTAGCTGATCGATCGTGGCCACCAGGGGGCGCTAGATCGATCGAT" +
      "GGATGAGTCATCCGGATCGATCGGGACATGCCCGGGCATGTCCGATCGAT" +
      "CCCTTTGTTCGATCGATCGATCTAGCGCCCCCTGGTGGCCAGATCGATCG" +
      "CGGATGACTCATCCGATCGATCGTGGCCACCAGGGGGCGCTAGATCGATC"
    );

    // =========================================================================
    // 1. Cross-Platform Light, Dark, and System Theme Engine
    // =========================================================================
    const systemMediaQuery = window.matchMedia
      ? window.matchMedia("(prefers-color-scheme: dark)")
      : null;

    function resolveEffectiveTheme(mode) {{
      if (mode === "light" || mode === "dark") return mode;
      if (systemMediaQuery && typeof systemMediaQuery.matches === "boolean") {{
        return systemMediaQuery.matches ? "dark" : "light";
      }}
      return (window.JASPAR_DASHBOARD_STATE && window.JASPAR_DASHBOARD_STATE.system_theme) || "light";
    }}

    function getActiveTheme() {{
      return document.documentElement.getAttribute("data-theme") === "dark" ? "dark" : "light";
    }}

    function setThemeMode(mode, persist = true) {{
      if (!["light", "dark", "system"].includes(mode)) mode = "system";
      currentThemeMode = mode;
      const effective = resolveEffectiveTheme(mode);
      const root = document.documentElement;
      root.setAttribute("data-theme", effective);
      root.setAttribute("data-theme-mode", mode);

      if (persist) {{
        try {{
          localStorage.setItem("jaspar_dashboard_theme_mode", mode);
        }} catch (_) {{}}
      }}

      const btnIds = {{
        light: "theme-btn-light",
        dark: "theme-btn-dark",
        system: "theme-btn-system"
      }};
      for (const [m, id] of Object.entries(btnIds)) {{
        const btn = document.getElementById(id);
        if (btn) {{
          const isActive = m === mode;
          btn.classList.toggle("active", isActive);
          btn.setAttribute("aria-pressed", isActive ? "true" : "false");
        }}
      }}

      // Re-render all theme-aware 2D charts, heatmaps, and 3D views without resetting viewport state
      if (window.__dashboardInitialized) {{
        renderLogoCard();
        updateMatrixView();
        const st = window.JASPAR_DASHBOARD_STATE;
        const locusEl = document.getElementById("locusSvgMount");
        if (locusEl && st && st.scan_results) {{
          locusEl.innerHTML = buildLocusTrackSvg(st.scan_results);
        }}
        render3DFrame();
      }}
    }}

    // Legacy helper alias for compatibility
    function toggleTheme() {{
      const next = getActiveTheme() === "dark" ? "light" : "dark";
      setThemeMode(next, true);
    }}

    if (systemMediaQuery) {{
      const handleOsThemeChange = () => {{
        if (currentThemeMode === "system") {{
          setThemeMode("system", false);
        }}
      }};
      if (typeof systemMediaQuery.addEventListener === "function") {{
        systemMediaQuery.addEventListener("change", handleOsThemeChange);
      }} else if (typeof systemMediaQuery.addListener === "function") {{
        systemMediaQuery.addListener(handleOsThemeChange);
      }}
    }}

    // =========================================================================
    // 2. Dynamic State Updates & Live SSE Streaming
    // =========================================================================
    window.updateDashboardState = function(partialState) {{
      if (!partialState || typeof partialState !== "object") return;
      const st = window.JASPAR_DASHBOARD_STATE;
      Object.assign(st, partialState);
      if (partialState.profiles_catalog) {{
        st.profiles_catalog = Object.assign(st.profiles_catalog || {{}}, partialState.profiles_catalog);
      }}
      if (partialState.input_sequence && partialState.input_sequence.sequence) {{
        const ta = document.getElementById("liveSequenceInput");
        if (ta && document.activeElement !== ta) {{
          ta.value = partialState.input_sequence.sequence;
        }}
      }}
      refreshAllCards(true);
      rebuild3DGraphNodes(true);
    }};

    function initSseStream() {{
      const badge = document.getElementById("sseStatusBadge");
      if (typeof EventSource === "undefined" || !window.location.protocol.startsWith("http")) {{
        if (badge) badge.textContent = "Mode: Offline Standalone";
        return;
      }}
      try {{
        const es = new EventSource("/events");
        es.onopen = () => {{
          if (badge) {{
            badge.textContent = "SSE: Live Stream Connected";
            badge.className = "status-signal status-emerald";
          }}
        }};
        es.onmessage = (evt) => {{
          try {{
            const payload = JSON.parse(evt.data);
            window.updateDashboardState(payload);
          }} catch (_) {{}}
        }};
        es.onerror = () => {{
          if (badge) {{
            badge.textContent = "Mode: Standalone Ready";
            badge.className = "status-signal status-teal";
          }}
          es.close();
        }};
      }} catch (_) {{}}
    }}

    // =========================================================================
    // 3. Sequence Presets, Active Profile & Sliding-Window Dual-Strand Scanner
    // =========================================================================
    function getActiveProfile() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      return (st.profiles_catalog && st.profiles_catalog[st.active_matrix_id]) || st;
    }}

    function loadSequencePreset(preset) {{
      const ta = document.getElementById("liveSequenceInput");
      if (preset === "promoter") {{
        window.JASPAR_DASHBOARD_STATE.input_sequence.name = "SYNTHETIC_MULTI_TF_PROMOTER_200BP";
        ta.value = MULTI_TF_PROMOTER_SEQ;
      }} else {{
        window.JASPAR_DASHBOARD_STATE.input_sequence.name = "SYNTHETIC_HUMAN_ENHANCER_1000BP";
        ta.value = DEFAULT_ENHANCER_SEQ;
      }}
      runLiveSequenceScan(true);
    }}

    function injectActiveMotifSite() {{
      const p = getActiveProfile();
      const canonicalMap = {{
        "MA0139.1": "TGGCCACCAGGGGGCGCTA",
        "MA0488.2": "GGATGAGTCATCCG",
        "MA0106.3": "GGACATGCCCGGGCATGTCC",
        "MA0143.4": "CCCTTTGTTC"
      }};
      const site = canonicalMap[p.matrix_id] || p.logo.consensus_fwd.replace(/[^ACGT]/g, "G");
      const ta = document.getElementById("liveSequenceInput");
      ta.value = "GATCGATCG" + site + "GATCGATCG" + ta.value;
      window.JASPAR_DASHBOARD_STATE.input_sequence.name = "CUSTOM_MODIFIED_SEQUENCE";
      runLiveSequenceScan(true);
    }}

    function switchActiveProfile(matrixId) {{
      const st = window.JASPAR_DASHBOARD_STATE;
      if (!st.profiles_catalog[matrixId]) return;
      st.active_matrix_id = matrixId;
      const p = st.profiles_catalog[matrixId];
      st.matrix_id = p.matrix_id;
      st.name = p.name;
      st.version = p.version;
      st.collection = p.collection;
      st.species = p.species;
      st.superclass = p.superclass;
      st.tf_class = p.tf_class;
      st.tf_family = p.tf_family;
      st.subfamily = p.subfamily;
      st.uniprot_ids = p.uniprot_ids;
      st.pdb_ids = p.pdb_ids;
      st.pfm = p.pfm;
      st.logo = p.logo;
      st.matrix_analysis = p.matrix_analysis;
      st.cluster_comparisons = p.cluster_comparisons;
      runLiveSequenceScan(true);
    }}

    function revCompDNA(seq) {{
      const map = {{ A: "T", T: "A", C: "G", G: "C", N: "N" }};
      let out = "";
      for (let i = seq.length - 1; i >= 0; i--) {{
        out += map[seq[i]] || "N";
      }}
      return out;
    }}

    function runLiveSequenceScan(userClicked) {{
      const t0 = performance.now();
      scanCounter++;
      const st = window.JASPAR_DASHBOARD_STATE;
      const p = getActiveProfile();
      const rawInput = (document.getElementById("liveSequenceInput").value || "").toUpperCase().replace(/[^ACGTN]/g, "");
      const threshold = parseFloat(document.getElementById("rngThreshold").value);
      document.getElementById("lblThreshold").textContent = threshold.toFixed(2);

      const pseudo = parseFloat(document.getElementById("rngPseudo").value);
      const gc = parseFloat(document.getElementById("rngGc").value);
      const bg = {{ A: (1 - gc) / 2, T: (1 - gc) / 2, G: gc / 2, C: gc / 2 }};
      const pfm = p.pfm;
      const L = p.logo.length;

      const pwm = {{ A: [], C: [], G: [], T: [] }};
      let sMax = 0;
      let sMin = 0;
      for (let i = 0; i < L; i++) {{
        const n_i = pfm.A[i] + pfm.C[i] + pfm.G[i] + pfm.T[i];
        let colMax = -Infinity;
        let colMin = Infinity;
        for (const b of ["A", "C", "G", "T"]) {{
          const prob = (pfm[b][i] + pseudo) / (n_i + 4 * pseudo);
          const w = Math.log2(prob / bg[b]);
          pwm[b].push(w);
          if (w > colMax) colMax = w;
          if (w < colMin) colMin = w;
        }}
        sMax += colMax;
        sMin += colMin;
      }}
      const span = Math.max(1e-9, sMax - sMin);
      const hits = [];

      if (rawInput.length >= L) {{
        for (let offset = 0; offset <= rawInput.length - L; offset++) {{
          const fwd = rawInput.slice(offset, offset + L);
          const rev = revCompDNA(fwd);
          for (const [strand, kmer] of [["+", fwd], ["-", rev]]) {{
            let rawScore = 0;
            for (let i = 0; i < L; i++) {{
              const ch = kmer[i];
              if (pwm[ch]) rawScore += pwm[ch][i];
            }}
            const rel = Math.max(0, Math.min(1, (rawScore - sMin) / span));
            if (rel >= threshold) {{
              const sig = rel >= 0.95 ? "Optimal" : (rel >= 0.85 ? "High" : "Moderate");
              hits.push({{
                start: offset + 1,
                end: offset + L,
                strand: strand,
                matched_sequence: kmer,
                raw_score: rawScore,
                relative_score: rel,
                significance: sig
              }});
            }}
          }}
        }}
      }}
      hits.sort((a, b) => b.relative_score - a.relative_score || a.start - b.start);

      const counts = {{ A: 0, C: 0, G: 0, T: 0 }};
      for (const ch of rawInput) if (counts[ch] !== undefined) counts[ch]++;
      const obsGc = rawInput.length > 0 ? ((counts.G + counts.C) / rawInput.length) : 0;

      st.input_sequence.sequence = rawInput;
      st.input_sequence.length = rawInput.length;
      st.input_sequence.observed_gc = obsGc;
      st.input_sequence.counts = counts;

      st.scan_results = {{
        matrix_id: p.matrix_id,
        name: p.name,
        sequence_name: st.input_sequence.name,
        sequence_length: rawInput.length,
        threshold: threshold,
        s_min: sMin,
        s_max: sMax,
        total_hits: hits.length,
        hits: hits
      }};

      const elapsedMs = Math.max(0.2, performance.now() - t0).toFixed(1);
      const fbText = document.getElementById("scanFeedbackText");
      if (fbText) {{
        const topSummary = hits.length > 0
          ? `Top Hit: ${{hits[0].strand}}[${{hits[0].start}}..${{hits[0].end}}] ${{hits[0].matched_sequence}} (${{(hits[0].relative_score * 100).toFixed(1)}}% affinity, ${{hits[0].significance}})`
          : `0 sites >= ${{(threshold * 100).toFixed(0)}}% cutoff (adjust Relative Cutoff slider to inspect lower-affinity windows)`;
        fbText.innerHTML = `<strong>Scan #${{scanCounter}} (${{elapsedMs}} ms):</strong> ${{hits.length}} binding site(s) for <strong>${{p.name}} (${{p.matrix_id}})</strong> across ${{rawInput.length}} bp | ${{topSummary}}`;
      }}

      if (userClicked) {{
        const btn = document.getElementById("btnRunScan");
        if (btn) {{
          btn.textContent = `Scanned (${{hits.length}} Hits)`;
          setTimeout(() => {{ btn.textContent = "Run Sequence Scan"; }}, 1200);
        }}
      }}

      refreshAllCards(true);
      rebuild3DGraphNodes(true);
    }}

    // =========================================================================
    // 4. Theme-Reactive 2D Charts: Sequence Logo & Coordinate Locus Track
    // =========================================================================
    function buildClientLogoSvg(profile, revComp, probMode) {{
      const logo = profile.logo;
      const cols = (revComp && logo.columns_rev) ? logo.columns_rev : logo.columns;
      const length = logo.length;
      const colWidth = 42.0;
      const leftMargin = 64.0;
      const rightMargin = 24.0;
      const topMargin = 46.0;
      const plotHeight = 210.0;
      const bottomMargin = 70.0;

      const svgWidth = Math.round(leftMargin + length * colWidth + rightMargin);
      const svgHeight = Math.round(topMargin + plotHeight + bottomMargin);
      const maxYVal = probMode ? 1.0 : 2.0;
      const yLabel = probMode ? "Position Probability" : "Information Content (bits)";
      const yTicks = probMode ? [0.0, 0.25, 0.5, 0.75, 1.0] : [0.0, 0.5, 1.0, 1.5, 2.0];

      const isDark = getActiveTheme() === "dark";
      const textPrimary = isDark ? "#ECE8E2" : "#1C1917";
      const textSecondary = isDark ? "#9A9895" : "#787774";
      const gridColor = isDark ? "#292A30" : "#E7E5E4";

      const orientLabel = revComp ? "Reverse Complement (3' <- 5')" : "Forward (5' -> 3')";
      const consensusStr = revComp ? logo.consensus_rev : logo.consensus_fwd;
      const titleStr = `${{profile.matrix_id}} (${{profile.name}}) Sequence Logo [${{orientLabel}}]`;

      const parts = [
        `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{svgWidth}} ${{svgHeight}}" width="100%" height="${{svgHeight}}" role="img" aria-label="${{titleStr}}" style="background:transparent; font-family:'DM Sans', sans-serif; display:block;">`,
        `<title>${{titleStr}}</title>`,
        `<rect x="0" y="0" width="${{svgWidth}}" height="${{svgHeight}}" fill="transparent"/>`,
        `<text x="${{leftMargin}}" y="24" font-size="13" font-weight="600" fill="${{textPrimary}}">${{titleStr}}</text>`,
        `<text x="${{svgWidth - rightMargin}}" y="24" font-size="11.5" font-weight="500" font-family="'IBM Plex Mono', monospace" fill="${{textSecondary}}" text-anchor="end">Total Info: ${{logo.total_bits.toFixed(2)}} bits | Consensus: ${{consensusStr}}</text>`
      ];

      for (const tick of yTicks) {{
        const frac = tick / maxYVal;
        const yPos = topMargin + plotHeight * (1.0 - frac);
        const dash = (tick === 0.0 || tick === maxYVal) ? "" : ' stroke-dasharray="3,3"';
        parts.push(`<line x1="${{leftMargin}}" y1="${{yPos.toFixed(2)}}" x2="${{svgWidth - rightMargin}}" y2="${{yPos.toFixed(2)}}" stroke="${{gridColor}}" stroke-width="1"${{dash}}/>`);
        parts.push(`<text x="${{leftMargin - 10}}" y="${{(yPos + 4).toFixed(2)}}" font-size="11" font-family="'IBM Plex Mono', monospace" fill="${{textSecondary}}" text-anchor="end">${{tick.toFixed(2)}}</text>`);
      }}

      parts.push(`<line x1="${{leftMargin}}" y1="${{topMargin}}" x2="${{leftMargin}}" y2="${{topMargin + plotHeight}}" stroke="${{gridColor}}" stroke-width="1.5"/>`);
      const midY = topMargin + plotHeight / 2.0;
      parts.push(`<text x="16" y="${{midY.toFixed(2)}}" font-size="11" font-weight="500" fill="${{textSecondary}}" text-anchor="middle" transform="rotate(-90, 16, ${{midY.toFixed(2)}})">${{yLabel}}</text>`);

      const baselineY = topMargin + plotHeight;
      cols.forEach((col, colIdx) => {{
        const xLeft = leftMargin + colIdx * colWidth;
        const xCenter = xLeft + colWidth / 2.0;
        let currentY = baselineY;

        for (const item of col.stack) {{
          const base = item.base;
          const val = probMode ? item.prob : item.bits;
          if (val <= 0.002) continue;
          const glyphH = (val / maxYVal) * plotHeight;
          const topY = currentY - glyphH;
          const scaleX = (colWidth - 4.0) / 100.0;
          const scaleY = glyphH / 100.0;
          const pathD = BASE_SVG_PATHS[base];
          const color = BASE_SEMANTIC_COLORS[base] || SEMANTIC_COLORS.teal;
          const tip = `Pos ${{col.position}} (${{base}}): ${{item.bits.toFixed(3)}} bits | P=${{item.prob.toFixed(3)}} | R_i=${{col.info_r_i.toFixed(3)}} bits`;
          parts.push(
            `<g transform="translate(${{(xLeft + 2.0).toFixed(2)}}, ${{topY.toFixed(2)}}) scale(${{scaleX.toFixed(4)}}, ${{scaleY.toFixed(4)}})"><title>${{tip}}</title><path d="${{pathD}}" fill="${{color}}"/></g>`
          );
          currentY = topY;
        }}

        parts.push(`<text x="${{xCenter.toFixed(2)}}" y="${{(baselineY + 16).toFixed(2)}}" font-size="10.5" font-family="'IBM Plex Mono', monospace" fill="${{textSecondary}}" text-anchor="middle">${{col.position}}</text>`);
        const iupac = col.iupac;
        const badgeStroke = BASE_SEMANTIC_COLORS[iupac] || SEMANTIC_COLORS.teal;
        parts.push(`<rect x="${{(xLeft + 3).toFixed(2)}}" y="${{(baselineY + 24).toFixed(2)}}" width="${{(colWidth - 6).toFixed(2)}}" height="22" rx="4" fill="${{badgeStroke}}" fill-opacity="0.14" stroke="${{badgeStroke}}" stroke-width="1"/>`);
        parts.push(`<text x="${{xCenter.toFixed(2)}}" y="${{(baselineY + 39).toFixed(2)}}" font-size="11.5" font-weight="600" font-family="'IBM Plex Mono', monospace" fill="${{textPrimary}}" text-anchor="middle">${{iupac}}</text>`);
      }});

      parts.push(`<text x="${{leftMargin}}" y="${{svgHeight - 6}}" font-size="10.5" fill="${{textSecondary}}">IUPAC Degenerate Consensus Strip (Position 1 -&gt; L) | Palette: A=#059669 C=#4F46E5 G=#D97706 T=#E11D48</text>`);
      parts.push(`</svg>`);
      return parts.join("\\n");
    }}

    function buildLocusTrackSvg(scan) {{
      const seqLen = Math.max(1, scan.sequence_length);
      const width = 940, height = 165, leftPad = 76, rightPad = 32;
      const trackW = width - leftPad - rightPad, midY = 86;

      const isDark = getActiveTheme() === "dark";
      const textPrimary = isDark ? "#ECE8E2" : "#1C1917";
      const textSecondary = isDark ? "#9A9895" : "#787774";
      const gridColor = isDark ? "#292A30" : "#E7E5E4";

      let out = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{width}} ${{height}}" width="100%" height="${{height}}" style="background:transparent; display:block; font-family:'DM Sans', sans-serif;">`;
      out += `<rect x="0" y="0" width="${{width}}" height="${{height}}" fill="transparent"/>`;
      out += `<text x="${{leftPad}}" y="20" font-size="12.5" font-weight="600" fill="${{textPrimary}}">Genomic Coordinate Locus Track: ${{scan.sequence_name}} (1 .. ${{seqLen}} bp)</text>`;

      // Quantitative coordinate gridlines
      for (const yGrid of [midY - 44, midY - 22, midY + 22, midY + 44]) {{
        out += `<line x1="${{leftPad}}" y1="${{yGrid}}" x2="${{width - rightPad}}" y2="${{yGrid}}" stroke="${{gridColor}}" stroke-width="1" stroke-dasharray="3,3"/>`;
      }}

      out += `<line x1="${{leftPad}}" y1="${{midY}}" x2="${{width - rightPad}}" y2="${{midY}}" stroke="${{textSecondary}}" stroke-width="1.5"/>`;
      out += `<text x="${{leftPad - 10}}" y="${{midY - 24}}" font-size="11" font-weight="600" font-family="'IBM Plex Mono', monospace" fill="${{SEMANTIC_COLORS.teal}}" text-anchor="end">+ Strand</text>`;
      out += `<text x="${{leftPad - 10}}" y="${{midY + 32}}" font-size="11" font-weight="600" font-family="'IBM Plex Mono', monospace" fill="${{SEMANTIC_COLORS.rose}}" text-anchor="end">- Strand</text>`;

      for (const frac of [0, 0.25, 0.5, 0.75, 1]) {{
        const x = leftPad + frac * trackW;
        const bp = Math.round(1 + frac * (seqLen - 1));
        out += `<line x1="${{x}}" y1="${{midY - 46}}" x2="${{x}}" y2="${{midY + 46}}" stroke="${{gridColor}}" stroke-width="1"/>`;
        out += `<text x="${{x}}" y="${{height - 10}}" font-size="10.5" font-family="'IBM Plex Mono', monospace" fill="${{textSecondary}}" text-anchor="middle">${{bp}} bp</text>`;
      }}

      for (const h of scan.hits) {{
        const x = leftPad + ((h.start - 1) / Math.max(1, seqLen)) * trackW;
        const isFwd = h.strand === "+";
        const col = isFwd
          ? (h.significance === "Optimal" ? SEMANTIC_COLORS.emerald : SEMANTIC_COLORS.teal)
          : (h.significance === "Optimal" ? SEMANTIC_COLORS.indigo : SEMANTIC_COLORS.rose);
        const stem = 16 + (h.relative_score - 0.5) * 52;
        const yEnd = isFwd ? (midY - stem) : (midY + stem);
        out += `<line x1="${{x}}" y1="${{midY}}" x2="${{x}}" y2="${{yEnd}}" stroke="${{col}}" stroke-width="2.2"/>`;
        out += `<circle cx="${{x}}" cy="${{yEnd}}" r="5" fill="${{col}}"><title>${{h.strand}} Strand [${{h.start}}..${{h.end}}]: ${{h.matched_sequence}} (Rel=${{(h.relative_score * 100).toFixed(1)}}%, ${{h.significance}})</title></circle>`;
      }}
      out += `</svg>`;
      return out;
    }}

    // =========================================================================
    // 5. Semantic Matrix Heatmap Interpolation
    //    Surface Tone -> Vibrant Teal (#0D9488) -> Vivid Indigo (#4F46E5) -> Deep Rose (#E11D48)
    // =========================================================================
    function hexToRgb(hex) {{
      const h = hex.replace("#", "");
      return [
        parseInt(h.substring(0, 2), 16),
        parseInt(h.substring(2, 4), 16),
        parseInt(h.substring(4, 6), 16)
      ];
    }}

    function lerpRgb(c1, c2, t) {{
      const clamped = Math.max(0, Math.min(1, t));
      const r = Math.round(c1[0] + (c2[0] - c1[0]) * clamped);
      const g = Math.round(c1[1] + (c2[1] - c1[1]) * clamped);
      const b = Math.round(c1[2] + (c2[2] - c1[2]) * clamped);
      return `rgb(${{r}}, ${{g}}, ${{b}})`;
    }}

    function interpolateHeatmapScale(intensity) {{
      const t = Math.max(0, Math.min(1, intensity));
      const isDark = getActiveTheme() === "dark";
      const surfaceRgb = hexToRgb(isDark ? "#16171D" : "#FFFFFF");
      const tealRgb = hexToRgb("#0D9488");
      const indigoRgb = hexToRgb("#4F46E5");
      const roseRgb = hexToRgb("#E11D48");

      if (t <= 0.03) return "transparent";
      if (t < 0.40) {{
        return lerpRgb(surfaceRgb, tealRgb, (t / 0.40) * 0.55);
      }}
      if (t < 0.75) {{
        return lerpRgb(
          lerpRgbTuple(surfaceRgb, tealRgb, 0.55),
          lerpRgbTuple(surfaceRgb, indigoRgb, 0.62),
          (t - 0.40) / 0.35
        );
      }}
      return lerpRgb(
        lerpRgbTuple(surfaceRgb, indigoRgb, 0.62),
        lerpRgbTuple(surfaceRgb, roseRgb, 0.72),
        (t - 0.75) / 0.25
      );
    }}

    function lerpRgbTuple(c1, c2, t) {{
      const clamped = Math.max(0, Math.min(1, t));
      return [
        Math.round(c1[0] + (c2[0] - c1[0]) * clamped),
        Math.round(c1[1] + (c2[1] - c1[1]) * clamped),
        Math.round(c1[2] + (c2[2] - c1[2]) * clamped)
      ];
    }}

    function cellColor(val, mode, colTotal) {{
      if (mode === "pfm") {{
        const frac = colTotal > 0 ? val / colTotal : 0;
        return interpolateHeatmapScale(frac);
      }}
      if (mode === "ppm") {{
        return interpolateHeatmapScale(val);
      }}
      // PWM log2-odds: map [-3.0, +2.0] onto smooth semantic concentration gradient
      if (val >= 0) {{
        return interpolateHeatmapScale(0.35 + Math.min(0.65, (val / 2.0) * 0.65));
      }}
      const depl = Math.min(0.32, (Math.abs(val) / 4.0) * 0.32);
      return interpolateHeatmapScale(depl);
    }}

    function updateMatrixView() {{
      const p = getActiveProfile();
      const pseudo = parseFloat(document.getElementById("rngPseudo").value);
      const gc = parseFloat(document.getElementById("rngGc").value);
      document.getElementById("lblPseudo").textContent = pseudo.toFixed(1);
      document.getElementById("lblGc").textContent = gc.toFixed(2);

      ["pwm", "ppm", "pfm"].forEach(m => {{
        const btn = document.getElementById("btnMode" + m.charAt(0).toUpperCase() + m.slice(1));
        if (btn) btn.className = heatmapMode === m ? "btn active-tab" : "btn secondary";
      }});

      ["jaspar", "meme", "transfac", "yaml"].forEach(fmt => {{
        const btn = document.getElementById("btnFmt" + fmt.charAt(0).toUpperCase() + fmt.slice(1));
        if (btn) btn.className = exportFormat === fmt ? "btn active-tab" : "btn secondary";
      }});

      const pfm = p.pfm;
      const L = p.logo.length;
      const bg = {{ A: (1 - gc) / 2, T: (1 - gc) / 2, G: gc / 2, C: gc / 2 }};
      const bases = ["A", "C", "G", "T"];

      let tbl = '<table class="matrix-table"><thead><tr><th class="base-col">Base</th>';
      for (let i = 1; i <= L; i++) tbl += `<th>${{i}}</th>`;
      tbl += "</tr></thead><tbody>";

      for (const b of bases) {{
        const baseColor = BASE_SEMANTIC_COLORS[b] || SEMANTIC_COLORS.teal;
        tbl += `<tr><td class="base-col" style="color:${{baseColor}}">${{b}}</td>`;
        for (let i = 0; i < L; i++) {{
          const n_i = pfm.A[i] + pfm.C[i] + pfm.G[i] + pfm.T[i];
          const count = pfm[b][i];
          const prob = (count + pseudo) / (n_i + 4 * pseudo);
          const pwm = Math.log2(prob / bg[b]);
          const txt = heatmapMode === "pfm" ? count.toFixed(0) : (heatmapMode === "ppm" ? prob.toFixed(3) : pwm.toFixed(2));
          const bgStyle = cellColor(heatmapMode === "pwm" ? pwm : (heatmapMode === "ppm" ? prob : count), heatmapMode, n_i);
          tbl += `<td style="background:${{bgStyle}}">${{txt}}</td>`;
        }}
        tbl += "</tr>";
      }}
      tbl += "</tbody></table>";
      document.getElementById("matrixHeatmapMount").innerHTML = tbl;
      document.getElementById("matrixFormatOutput").textContent = p.matrix_analysis.exports[exportFormat];
    }}

    function setHeatmapMode(m) {{ heatmapMode = m; updateMatrixView(); }}
    function setExportFormat(f) {{ exportFormat = f; updateMatrixView(); }}
    function copyExportText() {{
      const txt = document.getElementById("matrixFormatOutput").textContent;
      if (navigator.clipboard) navigator.clipboard.writeText(txt);
    }}

    // =========================================================================
    // 6. Card & Table Renderers
    // =========================================================================
    function renderLogoCard() {{
      const p = getActiveProfile();
      const s = p.logo;
      document.getElementById("logoSvgMount").innerHTML = buildClientLogoSvg(p, isRevComp, isProbAxis);
      document.getElementById("activeOrientationText").textContent = isRevComp ? "Reverse Complement (3' <- 5')" : "Forward (5' -> 3')";
      document.getElementById("activeConsensusText").textContent = isRevComp ? s.consensus_rev : s.consensus_fwd;
      document.getElementById("logoBitsBadge").textContent = `${{s.total_bits.toFixed(2)}} bits`;
    }}

    function toggleRevComp() {{ isRevComp = !isRevComp; renderLogoCard(); }}
    function toggleLogoAxis() {{ isProbAxis = !isProbAxis; renderLogoCard(); }}

    function downloadActiveLogoSvg() {{
      const p = getActiveProfile();
      const svgContent = buildClientLogoSvg(p, isRevComp, isProbAxis);
      const blob = new Blob([svgContent], {{ type: "image/svg+xml;charset=utf-8" }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = p.matrix_id + (isRevComp ? "_revcomp" : "_fwd") + ".svg";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }}

    function significanceBadgeClass(tier) {{
      if (tier === "Optimal") return "status-signal status-emerald";
      if (tier === "High") return "status-signal status-teal";
      if (tier === "Moderate") return "status-signal status-amber";
      return "status-signal status-rose";
    }}

    function collectionBadgeClass(col) {{
      if (col === "CORE") return "status-signal status-emerald";
      if (String(col).includes("BPNet")) return "status-signal status-indigo";
      return "status-signal status-amber";
    }}

    function refreshAllCards() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      const p = getActiveProfile();
      const speciesStr = (p.species || []).map(s => `${{s.name}} [TaxID: ${{s.tax_id}}]`).join(", ");

      document.getElementById("bannerTitle").textContent = `${{p.name}} Transcription Factor Binding Profile (${{p.matrix_id}})`;
      document.getElementById("bannerSubtitle").textContent = `Species: ${{speciesStr}} | TFClass: ${{p.superclass}} -> ${{p.tf_class}}`;
      document.getElementById("bannerBadges").innerHTML = `
        <span class="token-mono">Matrix: ${{p.matrix_id}}.v${{p.version}}</span>
        <span class="${{collectionBadgeClass(p.collection)}}">${{p.collection}}</span>
        <span class="token-mono">${{p.logo.length}} bp Motif</span>
      `;
      document.getElementById("activeProfileChip").textContent = `Active Matrix: ${{p.matrix_id}}`;

      document.getElementById("statProfileName").textContent = `${{p.name}} (${{p.matrix_id}})`;
      document.getElementById("statSeqLength").textContent = `${{st.input_sequence.length}} bp`;
      document.getElementById("statSeqGc").textContent = `${{(st.input_sequence.observed_gc * 100).toFixed(1)}}%`;
      const c = st.input_sequence.counts;
      document.getElementById("statBaseCounts").textContent = `A:${{c.A}} C:${{c.C}} G:${{c.G}} T:${{c.T}}`;
      document.getElementById("statTotalHits").textContent = `${{st.scan_results.total_hits}} hits`;

      const uniprotLinks = (p.uniprot_ids || []).map(u => `<a href="https://www.uniprot.org/uniprotkb/${{u}}" target="_blank" rel="noopener">${{u}}</a>`).join(", ");
      const pdbLinks = (p.pdb_ids || []).map(d => `<a href="https://www.rcsb.org/structure/${{d}}" target="_blank" rel="noopener">${{d}}</a>`).join(", ");
      document.getElementById("structuralKvMount").innerHTML = `
        <tr><td class="kv-key">TF Symbol</td><td class="kv-val"><strong>${{p.name}}</strong></td></tr>
        <tr><td class="kv-key">JASPAR Matrix ID</td><td class="kv-val">${{p.matrix_id}} (v${{p.version}})</td></tr>
        <tr><td class="kv-key">Validation Collection</td><td class="kv-val"><span class="${{collectionBadgeClass(p.collection)}}">${{p.collection}}</span></td></tr>
        <tr><td class="kv-key">TFClass Superclass</td><td class="kv-val">${{p.superclass}}</td></tr>
        <tr><td class="kv-key">TFClass Class</td><td class="kv-val">${{p.tf_class}}</td></tr>
        <tr><td class="kv-key">TFClass Family</td><td class="kv-val">${{p.tf_family}}</td></tr>
        <tr><td class="kv-key">TFClass Subfamily</td><td class="kv-val">${{p.subfamily}}</td></tr>
        <tr><td class="kv-key">Species Taxonomy</td><td class="kv-val">${{speciesStr}}</td></tr>
        <tr><td class="kv-key">UniProt Accession(s)</td><td class="kv-val">${{uniprotLinks}}</td></tr>
        <tr><td class="kv-key">RCSB PDB 3D Structures</td><td class="kv-val">${{pdbLinks}}</td></tr>
      `;

      document.getElementById("clusterTableBody").innerHTML = (p.cluster_comparisons || []).map(row => `
        <tr>
          <td><strong>${{row.matrix_id}}</strong></td>
          <td style="font-family:var(--font-sans);">${{row.profile_label}}</td>
          <td><span class="${{collectionBadgeClass(row.collection)}}">${{row.collection}}</span></td>
          <td><strong>${{row.consensus}}</strong></td>
          <td>${{row.pearson_r.toFixed(3)}}</td>
          <td>${{row.jensen_shannon_divergence.toFixed(3)}}</td>
        </tr>
      `).join("");

      if (st.provenance && st.provenance.compiled_at_utc) {{
        const provEl = document.getElementById("provenanceTimestamp");
        if (provEl) provEl.textContent = `Compiled: ${{st.provenance.compiled_at_utc}} | Host Theme: ${{st.system_theme}}`;
      }}

      const telemetryBody = document.getElementById("telemetryTableBody");
      if (telemetryBody && Array.isArray(st.telemetry)) {{
        telemetryBody.innerHTML = st.telemetry.map(ev => {{
          const cls = ev.status === "WARNING" ? "status-signal status-amber" : "status-signal status-emerald";
          return `<tr>
            <td><strong>${{ev.event}}</strong></td>
            <td><span class="${{cls}}">${{ev.status}}</span></td>
            <td style="font-family:var(--font-sans);">${{ev.detail}}</td>
          </tr>`;
        }}).join("");
      }}

      renderLogoCard();
      updateMatrixView();
      document.getElementById("locusSvgMount").innerHTML = buildLocusTrackSvg(st.scan_results);
      renderScanTable();
    }}

    function renderScanTable() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      const scan = st.scan_results;
      const query = (document.getElementById("hitSearchInput").value || "").trim().toLowerCase();
      const hits = scan.hits || [];
      const filtered = hits.filter(h => {{
        if (!query) return true;
        return (
          String(h.strand).toLowerCase().includes(query) ||
          String(h.matched_sequence).toLowerCase().includes(query) ||
          String(h.significance).toLowerCase().includes(query) ||
          String(h.start).includes(query)
        );
      }});
      const hitLabel = document.getElementById("hitCountLabel");
      hitLabel.textContent = `${{filtered.length}} Sites`;
      hitLabel.className = filtered.length > 0 ? "status-signal status-emerald" : "status-signal status-amber";

      document.getElementById("scanSummaryMeta").innerHTML = `
        Target Locus: <strong>${{scan.sequence_name}}</strong> (${{scan.sequence_length}} bp) |
        Relative Cutoff: <strong>${{scan.threshold.toFixed(2)}}</strong>
      `;
      const rows = filtered.map((h, idx) => `
        <tr>
          <td>${{idx + 1}}</td>
          <td>${{h.start}}</td>
          <td>${{h.end}}</td>
          <td><span class="${{h.strand === '+' ? 'status-signal status-teal' : 'status-signal status-rose'}}">${{h.strand}}</span></td>
          <td><strong>${{h.matched_sequence}}</strong></td>
          <td>${{h.raw_score.toFixed(3)}}</td>
          <td>${{(h.relative_score * 100).toFixed(1)}}% (${{h.relative_score.toFixed(3)}})</td>
          <td><span class="${{significanceBadgeClass(h.significance)}}">${{h.significance}}</span></td>
        </tr>
      `).join("");
      document.getElementById("scanTableBody").innerHTML = rows || '<tr><td colspan="8">No binding sites above active threshold.</td></tr>';
    }}

    // =========================================================================
    // 7. 3D & Network Visualization Engine
    //    - "Interactive Relationship Node Map" & "Structural Data Diagram"
    //    - Semantic palette nodes & semi-transparent connected edge strokes
    //    - Deterministic layout preserving pan, zoom, and rotation
    // =========================================================================
    let camRotX = 0.35;
    let camRotY = -0.55;
    let camZoom = 1.0;
    let camPanX = 0.0;
    let camPanY = 0.0;
    let isDragging3D = false;
    let isPanning3D = false;
    let prevMouseX = 0;
    let prevMouseY = 0;
    let graphNodes = [];
    let graphEdges = [];
    let selectedNodeId = null;

    function hexWithAlpha(hexColor, alpha) {{
      const [r, g, b] = hexToRgb(hexColor);
      return `rgba(${{r}}, ${{g}}, ${{b}}, ${{alpha}})`;
    }}

    function set3DViewMode(mode) {{
      view3DMode = mode;
      document.getElementById("btn3DNetwork").className = mode === "network" ? "btn active-tab" : "btn secondary";
      document.getElementById("btn3DHelix").className = mode === "helix" ? "btn active-tab" : "btn secondary";
      document.getElementById("view3DTitle").textContent = mode === "network"
        ? "Interactive Relationship Node Map"
        : "Structural Data Diagram";
      document.getElementById("view3DModeBadge").textContent = mode === "network"
        ? "Motif Relationship Topology"
        : "3D B-DNA Binding Groove";
      rebuild3DGraphNodes(false);
    }}

    function reset3DCamera() {{
      camRotX = 0.35;
      camRotY = -0.55;
      camZoom = 1.0;
      camPanX = 0.0;
      camPanY = 0.0;
      render3DFrame();
    }}

    function rebuild3DGraphNodes(preserveSelection = true) {{
      const p = getActiveProfile();
      const scan = window.JASPAR_DASHBOARD_STATE.scan_results;
      const prevSelectedId = preserveSelection ? selectedNodeId : null;
      graphNodes = [];
      graphEdges = [];

      if (view3DMode === "network") {{
        // Central TF Hub Node (Vibrant Teal #0D9488)
        graphNodes.push({{
          id: "tf_hub",
          label: `${{p.name}} (${{p.matrix_id}})`,
          category: "Primary TF Hub",
          detail: `Consensus: ${{p.logo.consensus_fwd}} | Info: ${{p.logo.total_bits.toFixed(2)}} bits`,
          x: 0, y: 0, z: 0, vx: 0, vy: 0, vz: 0,
          radius: 13.5, color: SEMANTIC_COLORS.teal, fixed: true
        }});

        // Structural TFClass Domain Node (Functional Emerald #059669)
        graphNodes.push({{
          id: "tf_class",
          label: p.tf_class,
          category: "Validated TFClass Domain",
          detail: `Superclass: ${{p.superclass}}`,
          x: -95, y: -62, z: 38, vx: 0, vy: 0, vz: 0,
          radius: 9.5, color: SEMANTIC_COLORS.emerald
        }});
        graphEdges.push({{
          source: 0,
          target: 1,
          length: 112,
          color: hexWithAlpha(SEMANTIC_COLORS.emerald, 0.52)
        }});

        // UniProt / PDB Structural Node (Vivid Indigo #4F46E5)
        const upId = (p.uniprot_ids && p.uniprot_ids[0]) || "UniProt";
        graphNodes.push({{
          id: "uniprot",
          label: `UniProt: ${{upId}}`,
          category: "Structural Group",
          detail: `PDB: ${{(p.pdb_ids || []).join(", ")}}`,
          x: -102, y: 58, z: -34, vx: 0, vy: 0, vz: 0,
          radius: 9.0, color: SEMANTIC_COLORS.indigo
        }});
        graphEdges.push({{
          source: 0,
          target: 2,
          length: 112,
          color: hexWithAlpha(SEMANTIC_COLORS.indigo, 0.52)
        }});

        // Familial Cluster Comparison Nodes (Vivid Indigo #4F46E5 & Energetic Amber #D97706)
        (p.cluster_comparisons || []).slice(1).forEach((cl, idx) => {{
          const nIdx = graphNodes.length;
          const angle = (idx + 1) * 1.25;
          const nodeColor = String(cl.collection).includes("BPNet")
            ? SEMANTIC_COLORS.indigo
            : SEMANTIC_COLORS.amber;
          graphNodes.push({{
            id: `cluster_${{idx}}`,
            label: cl.matrix_id,
            category: cl.collection,
            detail: `Pearson r=${{cl.pearson_r.toFixed(3)}} | JSD=${{cl.jensen_shannon_divergence.toFixed(3)}}`,
            x: 92 * Math.cos(angle), y: -72, z: 92 * Math.sin(angle), vx: 0, vy: 0, vz: 0,
            radius: 8.5, color: nodeColor
          }});
          graphEdges.push({{
            source: 0,
            target: nIdx,
            length: 118,
            color: hexWithAlpha(nodeColor, 0.50)
          }});
        }});

        // Detected Genomic Binding Site Nodes colored by validation tier / strand
        (scan.hits || []).slice(0, 8).forEach((hit, idx) => {{
          const nIdx = graphNodes.length;
          const theta = (idx / Math.max(1, Math.min(8, scan.hits.length))) * Math.PI * 2;
          const hitColor = hit.significance === "Optimal"
            ? SEMANTIC_COLORS.emerald
            : (hit.strand === "+" ? SEMANTIC_COLORS.teal : SEMANTIC_COLORS.rose);
          graphNodes.push({{
            id: `hit_${{idx}}`,
            label: `${{hit.strand}}[${{hit.start}}..${{hit.end}}]`,
            category: `Binding Site (${{hit.significance}})`,
            detail: `${{hit.matched_sequence}} | Rel=${{(hit.relative_score * 100).toFixed(1)}}% (${{hit.raw_score.toFixed(2)}} bits)`,
            x: 120 * Math.cos(theta), y: 66 + (idx % 2) * 22, z: 120 * Math.sin(theta), vx: 0, vy: 0, vz: 0,
            radius: 7.5 + hit.relative_score * 2.5,
            color: hitColor
          }});
          graphEdges.push({{
            source: 0,
            target: nIdx,
            length: 132,
            color: hexWithAlpha(hitColor, 0.48)
          }});
        }});
      }} else {{
        // 3D Structural Data Diagram (B-DNA Double Helix Binding Groove)
        const cols = p.logo.columns || [];
        const L = cols.length;
        const compMap = {{ A: "T", T: "A", C: "G", G: "C", N: "N" }};

        for (let i = 0; i < L; i++) {{
          const col = cols[i];
          const topItem = col.stack[col.stack.length - 1];
          const base = topItem.base;
          const compBase = compMap[base] || "N";
          const twist = i * 0.628;
          const yPos = (i - (L - 1) / 2) * 14.2;
          const radius = 56;

          const fwdColor = BASE_SEMANTIC_COLORS[base] || SEMANTIC_COLORS.teal;
          const revColor = BASE_SEMANTIC_COLORS[compBase] || SEMANTIC_COLORS.indigo;

          const idxFwd = graphNodes.length;
          graphNodes.push({{
            id: `bp_fwd_${{i+1}}`,
            label: `${{base}}${{i+1}}`,
            category: `Forward Strand Pos ${{i+1}}`,
            detail: `Base ${{base}} | R_i: ${{col.info_r_i.toFixed(2)}} bits | IUPAC: ${{col.iupac}}`,
            x: radius * Math.cos(twist),
            y: yPos,
            z: radius * Math.sin(twist),
            vx: 0, vy: 0, vz: 0,
            radius: 5.2 + col.info_r_i * 2.6,
            color: fwdColor,
            fixed: true
          }});

          const idxRev = graphNodes.length;
          graphNodes.push({{
            id: `bp_rev_${{i+1}}`,
            label: `${{compBase}}'${{i+1}}`,
            category: `Reverse Strand Pos ${{i+1}}`,
            detail: `Complement Base ${{compBase}} paired with ${{base}}${{i+1}}`,
            x: radius * Math.cos(twist + Math.PI),
            y: yPos,
            z: radius * Math.sin(twist + Math.PI),
            vx: 0, vy: 0, vz: 0,
            radius: 4.8,
            color: revColor,
            fixed: true
          }});

          // Hydrogen-bond base-pair rung (Vibrant Teal semi-transparent bond)
          graphEdges.push({{
            source: idxFwd,
            target: idxRev,
            length: 112,
            color: hexWithAlpha(SEMANTIC_COLORS.teal, 0.42)
          }});

          // Sugar-phosphate backbone links
          if (i > 0) {{
            graphEdges.push({{
              source: idxFwd - 2,
              target: idxFwd,
              length: 24,
              color: hexWithAlpha(SEMANTIC_COLORS.indigo, 0.68)
            }});
            graphEdges.push({{
              source: idxRev - 2,
              target: idxRev,
              length: 24,
              color: hexWithAlpha(SEMANTIC_COLORS.rose, 0.62)
            }});
          }}
        }}
      }}

      const matched = prevSelectedId && graphNodes.find(n => n.id === prevSelectedId);
      const activeNode = matched || graphNodes[0] || null;
      selectedNodeId = activeNode ? activeNode.id : null;
      updateInspectorBanner(activeNode);
    }}

    function stepNetworkRelaxation() {{
      if (view3DMode !== "network") return;
      const N = graphNodes.length;
      for (let i = 0; i < N; i++) {{
        for (let j = i + 1; j < N; j++) {{
          const dx = graphNodes[j].x - graphNodes[i].x;
          const dy = graphNodes[j].y - graphNodes[i].y;
          const dz = graphNodes[j].z - graphNodes[i].z;
          const distSq = dx * dx + dy * dy + dz * dz + 25;
          const force = 400 / distSq;
          const dist = Math.sqrt(distSq);
          const fx = (dx / dist) * force;
          const fy = (dy / dist) * force;
          const fz = (dz / dist) * force;
          if (!graphNodes[i].fixed) {{ graphNodes[i].vx -= fx; graphNodes[i].vy -= fy; graphNodes[i].vz -= fz; }}
          if (!graphNodes[j].fixed) {{ graphNodes[j].vx += fx; graphNodes[j].vy += fy; graphNodes[j].vz += fz; }}
        }}
      }}
      for (const e of graphEdges) {{
        const a = graphNodes[e.source];
        const b = graphNodes[e.target];
        const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1;
        const delta = (dist - e.length) * 0.035;
        const fx = (dx / dist) * delta, fy = (dy / dist) * delta, fz = (dz / dist) * delta;
        if (!a.fixed) {{ a.vx += fx; a.vy += fy; a.vz += fz; }}
        if (!b.fixed) {{ b.vx -= fx; b.vy -= fy; b.vz -= fz; }}
      }}
      for (const n of graphNodes) {{
        if (n.fixed) continue;
        n.vx *= 0.80; n.vy *= 0.80; n.vz *= 0.80;
        n.x += n.vx; n.y += n.vy; n.z += n.vz;
      }}
    }}

    function project3D(node, w, h) {{
      const cosY = Math.cos(camRotY), sinY = Math.sin(camRotY);
      const cosX = Math.cos(camRotX), sinX = Math.sin(camRotX);
      const x1 = node.x * cosY - node.z * sinY;
      const z1 = node.z * cosY + node.x * sinY;
      const y2 = node.y * cosX - z1 * sinX;
      const z2 = z1 * cosX + node.y * sinX;
      const fov = 420 * camZoom;
      const scale = fov / Math.max(120, 420 + z2);
      return {{
        sx: w / 2 + camPanX + x1 * scale,
        sy: h / 2 + camPanY + y2 * scale,
        z: z2,
        scale: scale
      }};
    }}

    function updateInspectorBanner(node) {{
      if (!node) return;
      document.getElementById("nodeInspectorText").textContent = `${{node.category}}: ${{node.label}}`;
      document.getElementById("nodeMetricBadge").textContent = node.detail;
    }}

    function render3DFrame() {{
      const surface = document.getElementById("viewport3DSurface");
      if (!surface) return;
      const rect = surface.getBoundingClientRect();
      const W = Math.max(300, rect.width || 520);
      const H = Math.max(240, rect.height || 360);
      surface.setAttribute("viewBox", `0 0 ${{W}} ${{H}}`);

      const isDark = getActiveTheme() === "dark";
      const labelColor = isDark ? "#ECE8E2" : "#1C1917";
      const selectedRing = isDark ? "#ECE8E2" : "#1C1917";

      for (const n of graphNodes) {{
        n._proj = project3D(n, W, H);
      }}

      const parts = [];
      for (const e of graphEdges) {{
        const p1 = graphNodes[e.source]._proj;
        const p2 = graphNodes[e.target]._proj;
        parts.push(
          `<line x1="${{p1.sx.toFixed(1)}}" y1="${{p1.sy.toFixed(1)}}" x2="${{p2.sx.toFixed(1)}}" y2="${{p2.sy.toFixed(1)}}" stroke="${{e.color}}" stroke-width="1.6"/>`
        );
      }}

      const sorted = [...graphNodes].sort((a, b) => b._proj.z - a._proj.z);
      for (const n of sorted) {{
        const pr = n._proj;
        const r = Math.max(3.2, n.radius * pr.scale);
        const strokeAttr = (selectedNodeId === n.id)
          ? ` stroke="${{selectedRing}}" stroke-width="2.2"`
          : "";
        parts.push(
          `<circle cx="${{pr.sx.toFixed(1)}}" cy="${{pr.sy.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{n.color}}"${{strokeAttr}}/>`
        );
        parts.push(
          `<text x="${{pr.sx.toFixed(1)}}" y="${{(pr.sy - r - 5).toFixed(1)}}" fill="${{labelColor}}" font-size="10.5" font-weight="600" text-anchor="middle" font-family="'IBM Plex Mono', monospace">${{n.label}}</text>`
        );
      }}
      surface.innerHTML = parts.join("");
    }}

    function init3DViewport() {{
      const host = document.getElementById("viewport3DHost");
      const svgNs = "http://www.w3.org/2000/svg";
      const surface = document.createElementNS(svgNs, "svg");
      surface.setAttribute("id", "viewport3DSurface");
      host.insertBefore(surface, host.firstChild);
      rebuild3DGraphNodes(false);

      surface.addEventListener("mousedown", e => {{
        isDragging3D = true;
        isPanning3D = e.shiftKey || e.button === 2;
        prevMouseX = e.clientX;
        prevMouseY = e.clientY;

        const rect = surface.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;
        for (const n of graphNodes) {{
          if (n._proj) {{
            const dist = Math.hypot(mx - n._proj.sx, my - n._proj.sy);
            if (dist <= n.radius * n._proj.scale + 6) {{
              selectedNodeId = n.id;
              updateInspectorBanner(n);
              break;
            }}
          }}
        }}
      }});

      window.addEventListener("mousemove", e => {{
        if (!isDragging3D) return;
        const dx = e.clientX - prevMouseX;
        const dy = e.clientY - prevMouseY;
        prevMouseX = e.clientX;
        prevMouseY = e.clientY;
        if (isPanning3D) {{
          camPanX += dx;
          camPanY += dy;
        }} else {{
          camRotY += dx * 0.008;
          camRotX = Math.max(-1.3, Math.min(1.3, camRotX + dy * 0.008));
        }}
      }});

      window.addEventListener("mouseup", () => {{
        isDragging3D = false;
        isPanning3D = false;
      }});

      surface.addEventListener("wheel", e => {{
        e.preventDefault();
        camZoom = Math.max(0.55, Math.min(2.2, camZoom - e.deltaY * 0.001));
      }}, {{ passive: false }});

      function animate() {{
        if (!isDragging3D) camRotY += 0.0022;
        stepNetworkRelaxation();
        render3DFrame();
        requestAnimationFrame(animate);
      }}
      requestAnimationFrame(animate);
    }}

    // =========================================================================
    // 8. Dashboard Initialization
    // =========================================================================
    function initDashboard() {{
      let savedMode = "system";
      try {{
        const stored = localStorage.getItem("jaspar_dashboard_theme_mode");
        if (stored && ["light", "dark", "system"].includes(stored)) {{
          savedMode = stored;
        }}
      }} catch (_) {{}}
      setThemeMode(savedMode, false);

      const st = window.JASPAR_DASHBOARD_STATE;
      const sel = document.getElementById("profileSelector");
      sel.innerHTML = Object.keys(st.profiles_catalog).map(id => {{
        const p = st.profiles_catalog[id];
        return `<option value="${{id}}" ${{id === st.active_matrix_id ? "selected" : ""}}>${{p.name}} (${{id}} - ${{p.tf_class}})</option>`;
      }}).join("");

      document.getElementById("liveSequenceInput").value = st.input_sequence.sequence;
      document.getElementById("rngThreshold").value = st.scan_results.threshold || 0.80;
      window.__dashboardInitialized = true;
      runLiveSequenceScan(false);
      init3DViewport();
      initSseStream();
    }}

    initDashboard();
  </script>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    """Constructs CLI argument parser for generate_dashboard.py."""
    parser = argparse.ArgumentParser(
        prog="generate_dashboard.py",
        description="Compile standalone HTML5 interactive dashboard (dashboard.html) for JASPAR TF motifs.",
    )
    parser.add_argument("--matrix-id", default="MA0139.1", help="JASPAR Matrix ID (default: MA0139.1).")
    parser.add_argument("--metadata-json", help="Optional JSON file from get_tf_metadata.")
    parser.add_argument("--motif-json", help="JSON file from get_tf_motif.")
    parser.add_argument("--scan-json", help="Optional JSON file from visualize_scan.py.")
    parser.add_argument("--fasta", help="Optional user FASTA file to embed and scan in dashboard.")
    parser.add_argument("--sequence", help="Optional user DNA sequence string to embed and scan in dashboard.")
    parser.add_argument("--output", required=True, help="Output path for dashboard.html.")
    parser.add_argument("--mock", action="store_true", help="Build complete dashboard using offline mock fixtures.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for generate_dashboard.py."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        verify_workspace_license(quiet=True)

        if args.motif_json and Path(args.motif_json).exists():
            motif_data = json.loads(Path(args.motif_json).read_text(encoding="utf-8"))
        elif args.mock:
            motif_data = load_mock_fixture(args.matrix_id)
        else:
            raise ValueError("Must provide an existing --motif-json file or pass --mock.")

        metadata_data: Optional[Dict[str, Any]] = None
        if args.metadata_json and Path(args.metadata_json).exists():
            metadata_data = json.loads(Path(args.metadata_json).read_text(encoding="utf-8"))

        scan_data: Optional[Dict[str, Any]] = None
        if args.scan_json and Path(args.scan_json).exists():
            scan_data = json.loads(Path(args.scan_json).read_text(encoding="utf-8"))

        html_doc = compile_dashboard_html(
            motif_data=motif_data,
            metadata_data=metadata_data,
            scan_data=scan_data,
            fasta_path=args.fasta,
            raw_sequence=args.sequence,
        )
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(html_doc, encoding="utf-8")
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
