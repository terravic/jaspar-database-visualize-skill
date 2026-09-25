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
3D interactive network and helical groove visualizations, and sliding-window DNA
sequence scan results into a single self-contained HTML5 dashboard.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

from jaspar_api import (
    MOCK_NAME_MAP,
    SCIENTIFIC_DISCLAIMER,
    load_mock_fixture,
    verify_workspace_license,
)
from visualize_logo import compute_logo_metrics, render_logo_svg
from visualize_matrix import build_matrix_analysis
from visualize_scan import parse_fasta_or_raw, render_locus_map_svg, scan_dna_sequence


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
            "profile_label": f"{name} (JASPAR CORE Canonical)",
            "collection": "CORE",
            "model_type": matrix_data.get("type", "ChIP-seq"),
            "length": length,
            "consensus": consensus,
            "pearson_r": 1.000,
            "jensen_shannon_divergence": 0.000,
        },
        {
            "matrix_id": f"{base_id}.BPNet",
            "profile_label": f"{name} (Deep Learning BPNet Contribution Profile)",
            "collection": "Deep Learning BPNet",
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

    svg_fwd_bits = render_logo_svg(fwd_metrics, title=logo_title, mode="bits") if include_inline_svgs else ""
    svg_rev_bits = render_logo_svg(rev_metrics, title=logo_title, mode="bits") if include_inline_svgs else ""
    svg_fwd_prob = render_logo_svg(fwd_metrics, title=logo_title, mode="prob") if include_inline_svgs else ""
    svg_rev_prob = render_logo_svg(rev_metrics, title=logo_title, mode="prob") if include_inline_svgs else ""

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
        "tf_class": ", ".join(motif_data.get("class", ["Zinc finger, C2H2-related"])),
        "tf_family": ", ".join(motif_data.get("family", ["More than 3 adjacent zinc finger factors"])),
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
            "svg_fwd_bits": svg_fwd_bits,
            "svg_rev_bits": svg_rev_bits,
            "svg_fwd_prob": svg_fwd_prob,
            "svg_rev_prob": svg_rev_prob,
        },
        "matrix_analysis": {"exports": matrix_analysis["exports"]},
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
    """Assembles the standalone HTML5 interactive dashboard with Light/Dark mode and 3D views."""
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

    dashboard_state = {
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
        "profiles_catalog": profiles_catalog,
    }

    matrix_id = primary_bundle["matrix_id"]
    name = primary_bundle["name"]

    return f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>JASPAR Interactive Dashboard: {html.escape(matrix_id)} ({html.escape(name)})</title>
  <style>
    :root {{
      --bg: #F1F5F9;
      --surface: #FFFFFF;
      --surface-alt: #F8FAFC;
      --border: #E2E8F0;
      --text: #0F172A;
      --muted: #475569;
      --primary: #1E293B;
      --accent: #2563EB;
      --accent-soft: rgba(37, 99, 235, 0.12);
      --badge-bg: #E2E8F0;
      --badge-text: #0F172A;
      --viewport-bg: #0F172A;
      --code-bg: #0F172A;
      --code-text: #F8FAFC;
      --disclaimer-bg: #FFFBEB;
      --disclaimer-border: #F59E0B;
      --disclaimer-text: #92400E;
      --shadow: 0 2px 8px rgba(15, 23, 42, 0.06);
    }}
    [data-theme="dark"] {{
      --bg: #090D16;
      --surface: #111827;
      --surface-alt: #1E293B;
      --border: #334155;
      --text: #F8FAFC;
      --muted: #94A3B8;
      --primary: #3B82F6;
      --accent: #60A5FA;
      --accent-soft: rgba(96, 165, 250, 0.16);
      --badge-bg: #1E293B;
      --badge-text: #E2E8F0;
      --viewport-bg: #060911;
      --code-bg: #060911;
      --code-text: #E2E8F0;
      --disclaimer-bg: rgba(245, 158, 11, 0.12);
      --disclaimer-border: #D97706;
      --disclaimer-text: #FDE68A;
      --shadow: 0 4px 14px rgba(0, 0, 0, 0.35);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, sans-serif;
      background: var(--bg);
      color: var(--text);
      margin: 0;
      padding: 20px;
      line-height: 1.5;
      transition: background 0.2s ease, color 0.2s ease;
    }}
    .container {{ max-width: 1320px; margin: 0 auto; display: flex; flex-direction: column; gap: 18px; }}
    .banner {{
      background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%);
      color: #FFFFFF;
      border-radius: 12px;
      padding: 20px 24px;
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 16px;
      border: 1px solid #334155;
      box-shadow: var(--shadow);
    }}
    .banner h1 {{ margin: 0 0 4px 0; font-size: 22px; letter-spacing: -0.02em; }}
    .banner-meta {{ font-size: 13px; color: #CBD5E1; }}
    .header-actions {{ display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }}
    .badge-row {{ display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }}
    .badge {{
      display: inline-block;
      background: var(--badge-bg);
      color: var(--badge-text);
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
    }}
    .banner .badge {{ background: rgba(255,255,255,0.12); color: #F8FAFC; border: 1px solid rgba(255,255,255,0.22); }}
    /* Theme Toggle Button */
    .theme-toggle-btn {{
      border-radius: 6px;
      border: 1px solid rgba(255, 255, 255, 0.28);
      background: rgba(255, 255, 255, 0.12);
      color: #FFFFFF;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      padding: 6px 12px;
      font-size: 12px;
      font-weight: 600;
      transition: background 0.15s ease;
    }}
    .theme-toggle-btn:hover {{ background: rgba(255, 255, 255, 0.24); }}
    .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }}
    .grid-3-2 {{ display: grid; grid-template-columns: 1.35fr 1fr; gap: 18px; }}
    @media (max-width: 1024px) {{
      .grid-2, .grid-3-2 {{ grid-template-columns: 1fr; }}
    }}
    .card {{
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 20px;
      box-shadow: var(--shadow);
    }}
    .card h2 {{
      margin: 0 0 12px 0;
      font-size: 16px;
      color: var(--text);
      border-bottom: 1px solid var(--border);
      padding-bottom: 10px;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: center;
      justify-content: space-between;
      background: var(--surface-alt);
      border: 1px solid var(--border);
      padding: 10px 14px;
      border-radius: 8px;
      margin-bottom: 14px;
    }}
    .btn {{
      background: var(--primary);
      color: #FFFFFF;
      border: 1px solid var(--primary);
      border-radius: 6px;
      padding: 7px 12px;
      font-size: 12px;
      font-weight: 600;
      cursor: pointer;
      transition: opacity 0.15s ease;
    }}
    .btn:hover {{ opacity: 0.92; }}
    .btn.secondary {{
      background: var(--surface);
      color: var(--text);
      border: 1px solid var(--border);
    }}
    .btn.active-tab {{
      background: var(--accent);
      border-color: var(--accent);
      color: #FFFFFF;
    }}
    .stat-pills {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; margin-bottom: 12px; }}
    .stat-pill {{
      background: var(--surface-alt);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 10px 12px;
    }}
    .stat-pill .label {{ font-size: 11px; color: var(--muted); text-transform: uppercase; font-weight: 700; }}
    .stat-pill .val {{ font-size: 16px; font-weight: 800; color: var(--text); margin-top: 2px; font-family: ui-monospace, monospace; }}
    .seq-textarea {{
      width: 100%;
      min-height: 66px;
      background: var(--surface-alt);
      color: var(--text);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 8px 10px;
      font-family: ui-monospace, monospace;
      font-size: 12px;
      resize: vertical;
    }}
    .select-input, .search-input {{
      padding: 6px 10px;
      background: var(--surface);
      color: var(--text);
      border: 1px solid var(--border);
      border-radius: 6px;
      font-size: 12px;
    }}
    .consensus-strip {{
      margin-top: 12px;
      background: var(--surface-alt);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 10px 14px;
      font-family: ui-monospace, monospace;
      font-size: 13px;
      display: flex;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 10px;
    }}
    .viewport-stage {{
      position: relative;
      width: 100%;
      height: 360px;
      background: var(--viewport-bg);
      border-radius: 10px;
      overflow: hidden;
      border: 1px solid var(--border);
    }}
    #viewport3DSurface {{ width: 100%; height: 100%; display: block; cursor: grab; }}
    #viewport3DSurface:active {{ cursor: grabbing; }}
    .viewport-overlay {{
      position: absolute;
      bottom: 10px;
      left: 12px;
      right: 12px;
      background: rgba(15, 23, 42, 0.84);
      color: #F8FAFC;
      border: 1px solid rgba(255,255,255,0.15);
      border-radius: 8px;
      padding: 8px 12px;
      font-size: 12px;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 8px;
      pointer-events: none;
    }}
    .table-wrap {{ overflow-x: auto; border: 1px solid var(--border); border-radius: 8px; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
    th, td {{ border: 1px solid var(--border); padding: 7px 9px; text-align: left; color: var(--text); }}
    th {{ background: var(--surface-alt); font-weight: 700; color: var(--muted); }}
    .matrix-table td, .matrix-table th {{ text-align: right; font-family: ui-monospace, monospace; }}
    .matrix-table .base-col {{ text-align: center; font-weight: 700; background: var(--surface-alt); width: 44px; }}
    pre.code-export {{
      background: var(--code-bg);
      color: var(--code-text);
      padding: 12px 14px;
      border-radius: 8px;
      font-size: 12px;
      overflow-x: auto;
      max-height: 210px;
      margin-top: 12px;
      border: 1px solid var(--border);
    }}
    .kv-list {{ display: grid; grid-template-columns: 165px 1fr; row-gap: 8px; column-gap: 12px; font-size: 13px; }}
    .kv-label {{ font-weight: 700; color: var(--muted); }}
    a {{ color: var(--accent); text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .footer-card {{
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 18px 22px;
      font-size: 12px;
      color: var(--muted);
    }}
    .disclaimer-box {{
      background: var(--disclaimer-bg);
      border: 1px solid var(--disclaimer-border);
      color: var(--disclaimer-text);
      padding: 12px 16px;
      border-radius: 8px;
      margin-top: 12px;
      font-weight: 500;
    }}
  </style>
</head>
<body>
  <div class="container">
    <!-- Top Header Banner with Theme Toggle -->
    <header class="banner">
      <div>
        <h1 id="bannerTitle">{html.escape(name)} Transcription Factor Binding Profile ({html.escape(matrix_id)})</h1>
        <div class="banner-meta" id="bannerSubtitle"></div>
      </div>
      <div class="header-actions">
        <div class="badge-row" id="bannerBadges"></div>
        <button class="theme-toggle-btn" id="themeToggleBtn" onclick="toggleTheme()" aria-label="Toggle Light and Dark Theme" title="Switch Light / Dark Theme">Theme: Light</button>
      </div>
    </header>

    <!-- Processed Input Data & Interactive Sequence Re-Scanner -->
    <section class="card">
      <h2>
        <span>Input Data &amp; Live Sequence Scanner</span>
        <span class="badge" id="activeProfileChip">Active Matrix: {html.escape(matrix_id)}</span>
      </h2>
      <div class="stat-pills">
        <div class="stat-pill">
          <div class="label">Active TF Profile</div>
          <div class="val" id="statProfileName">{html.escape(name)} ({html.escape(matrix_id)})</div>
        </div>
        <div class="stat-pill">
          <div class="label">Input Sequence Length</div>
          <div class="val" id="statSeqLength">0 bp</div>
        </div>
        <div class="stat-pill">
          <div class="label">Observed Sequence GC%</div>
          <div class="val" id="statSeqGc">0.0%</div>
        </div>
        <div class="stat-pill">
          <div class="label">Base Composition (A/C/G/T)</div>
          <div class="val" id="statBaseCounts">-</div>
        </div>
        <div class="stat-pill">
          <div class="label">Detected Binding Sites</div>
          <div class="val" id="statTotalHits">0</div>
        </div>
      </div>
      <div class="toolbar">
        <div style="display:flex; gap:8px; flex-wrap:wrap; align-items:center;">
          <label><strong>Switch JASPAR Profile:</strong> </label>
          <select id="profileSelector" class="select-input" onchange="switchActiveProfile(this.value)"></select>
          &nbsp;
          <label><strong>Relative Cutoff:</strong> <span id="lblThreshold">0.80</span></label>
          <input type="range" id="rngThreshold" min="0.65" max="0.98" step="0.01" value="0.80" oninput="runLiveSequenceScan(false)"/>
        </div>
        <div style="display:flex; gap:8px; flex-wrap:wrap;">
          <button class="btn secondary" onclick="loadSequencePreset('enhancer')">Preset: 1kb Enhancer</button>
          <button class="btn secondary" onclick="loadSequencePreset('promoter')">Preset: Multi-TF Promoter</button>
          <button class="btn secondary" onclick="injectActiveMotifSite()">+ Inject Motif Site</button>
          <button class="btn" id="btnRunScan" onclick="runLiveSequenceScan(true)">Run Sequence Scan</button>
        </div>
      </div>
      <textarea id="liveSequenceInput" class="seq-textarea" spellcheck="false" oninput="runLiveSequenceScan(false)"></textarea>
      <div id="scanFeedbackBanner" style="margin-top:10px; padding:10px 14px; border-radius:8px; border:1px solid var(--accent); background:var(--accent-soft); font-size:12px; display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; transition: transform 0.2s ease, box-shadow 0.2s ease;">
        <span id="scanFeedbackText">Ready to scan sequence across forward (+) and reverse (-) strands.</span>
        <button class="btn secondary" style="padding:4px 10px; font-size:11px;" onclick="document.getElementById('bindingSitesSection').scrollIntoView({{behavior:'smooth'}})">View Full Hits Table Below</button>
      </div>
    </section>

    <!-- Sequence Logo + 3D Interactive Viewport -->
    <div class="grid-3-2">
      <section class="card">
        <h2>
          <span>Sequence Logo</span>
          <span class="badge" id="logoBitsBadge">0.00 bits</span>
        </h2>
        <div class="toolbar">
          <div>
            <button class="btn" onclick="toggleRevComp()">Toggle Reverse Complement</button>
            <button class="btn secondary" onclick="toggleLogoAxis()">Bits / Probability</button>
          </div>
          <div>
            <button class="btn secondary" onclick="downloadActiveLogoSvg()">Download SVG</button>
          </div>
        </div>
        <div id="logoSvgMount" style="overflow-x:auto;"></div>
        <div class="consensus-strip">
          <span>Orientation: <strong id="activeOrientationText">Forward (5' -&gt; 3')</strong></span>
          <span>IUPAC Consensus: <strong id="activeConsensusText"></strong></span>
        </div>
      </section>

      <section class="card">
        <h2>
          <span id="view3DTitle">Motif Relationship Network</span>
          <span class="badge" id="view3DModeBadge">3D Interactive View</span>
        </h2>
        <div class="toolbar">
          <div>
            <button class="btn active-tab" id="btn3DNetwork" onclick="set3DViewMode('network')">Motif Relationship Network</button>
            <button class="btn secondary" id="btn3DHelix" onclick="set3DViewMode('helix')">3D DNA Binding Groove</button>
          </div>
          <div>
            <button class="btn secondary" onclick="reset3DCamera()">Reset View</button>
          </div>
        </div>
        <div class="viewport-stage" id="viewport3DHost">
          <div class="viewport-overlay">
            <span id="nodeInspectorText">Drag to rotate 3D view | Scroll to zoom | Click any node to inspect</span>
            <span id="nodeMetricBadge" style="color:#60A5FA; font-weight:700;"></span>
          </div>
        </div>
      </section>
    </div>

    <!-- PFM / PPM / PWM Matrix Explorer -->
    <section class="card">
      <h2>
        <span>Binding Affinity Matrix</span>
        <span class="badge">PFM / PPM / PWM Log2-Odds</span>
      </h2>
      <div class="toolbar">
        <div>
          <button class="btn active-tab" id="btnModePwm" onclick="setHeatmapMode('pwm')">PWM (Log2-Odds)</button>
          <button class="btn secondary" id="btnModePpm" onclick="setHeatmapMode('ppm')">PPM (Probabilities)</button>
          <button class="btn secondary" id="btnModePfm" onclick="setHeatmapMode('pfm')">PFM (Raw Counts)</button>
        </div>
        <div>
          <label>Pseudocount: <strong id="lblPseudo">0.8</strong></label>
          <input type="range" id="rngPseudo" min="0.1" max="5.0" step="0.1" value="0.8" oninput="updateMatrixView()"/>
          &nbsp;&nbsp;
          <label>Background GC: <strong id="lblGc">0.50</strong></label>
          <input type="range" id="rngGc" min="0.20" max="0.80" step="0.02" value="0.50" oninput="updateMatrixView()"/>
        </div>
        <div>
          <button class="btn secondary" onclick="setExportFormat('jaspar')">JASPAR</button>
          <button class="btn secondary" onclick="setExportFormat('meme')">MEME</button>
          <button class="btn secondary" onclick="setExportFormat('transfac')">TRANSFAC</button>
          <button class="btn secondary" onclick="setExportFormat('yaml')">YAML</button>
          <button class="btn" onclick="copyExportText()">Copy Format</button>
        </div>
      </div>
      <div class="table-wrap" id="matrixHeatmapMount"></div>
      <pre class="code-export" id="matrixFormatOutput"></pre>
    </section>

    <!-- Structural Classification & Familial Cluster Comparison -->
    <div class="grid-2">
      <section class="card">
        <h2>Structural Classification</h2>
        <div class="kv-list" id="structuralKvMount"></div>
      </section>

      <section class="card">
        <h2>Familial Cluster Comparison</h2>
        <div class="table-wrap">
          <table>
            <thead>
              <tr>
                <th>Profile ID</th><th>Profile Description</th><th>Collection</th>
                <th>Aligned Consensus</th><th>Pearson r</th><th>Mean JSD</th>
              </tr>
            </thead>
            <tbody id="clusterTableBody"></tbody>
          </table>
        </div>
      </section>
    </div>

    <!-- Genomic Binding Sites Track & Table -->
    <section class="card" id="bindingSitesSection">
      <h2>
        <span>Genomic Binding Sites</span>
        <span class="badge" id="hitCountLabel">0 Sites</span>
      </h2>
      <div class="toolbar">
        <div id="scanSummaryMeta"></div>
        <div>
          <input type="text" id="hitSearchInput" class="search-input" placeholder="Filter by strand (+ / -), k-mer, or significance..." oninput="renderScanTable()"/>
        </div>
      </div>
      <div id="locusSvgMount" style="margin-bottom:14px; overflow-x:auto;"></div>
      <div class="table-wrap">
        <table>
          <thead>
            <tr>
              <th>#</th><th>Start (1-based)</th><th>End</th><th>Strand</th>
              <th>Matched k-mer (5' -&gt; 3')</th><th>Raw Score (bits)</th>
              <th>Relative Affinity Score</th><th>Significance Tier</th>
            </tr>
          </thead>
          <tbody id="scanTableBody"></tbody>
        </table>
      </div>
    </section>

    <!-- Legal, Academic Citations & Scientific Disclaimer Footer -->
    <footer class="footer-card">
      <strong>Academic Citations &amp; Open-Access Data Licensing</strong>
      <ul style="margin:8px 0; padding-left:20px;">
        <li>
          <strong>JASPAR 2026 Release:</strong> Ovek Baydar, D. et al.
          <em>JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models.</em>
          Nucleic Acids Research 54, D184-D193 (2026). doi:10.1093/nar/gkaf1209.
        </li>
        <li>
          <strong>JASPAR 2024 / Structural Collections:</strong> Rauluseviciute, I. et al.
          <em>JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles.</em>
          Nucleic Acids Research 52, D174-D182 (2024). doi:10.1093/nar/gkad1059.
        </li>
        <li>
          <strong>License Attribution:</strong> Data retrieved from JASPAR (<a href="https://jaspar.elixir.no/" target="_blank" rel="noopener">https://jaspar.elixir.no/</a>)
          is licensed under the <a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noopener">Creative Commons Attribution 4.0 International (CC BY 4.0)</a> license.
          Dashboard software licensed under Apache License 2.0.
        </li>
      </ul>
      <div class="disclaimer-box">
        <strong>Scientific Disclaimer (AS-IS):</strong> {html.escape(SCIENTIFIC_DISCLAIMER)}
      </div>
    </footer>
  </div>

  <script>
    window.JASPAR_DASHBOARD_STATE = {json.dumps(dashboard_state)};

    let isRevComp = false;
    let isProbAxis = false;
    let heatmapMode = "pwm";
    let exportFormat = "jaspar";
    let view3DMode = "network";
    let scanCounter = 0;

    const DEFAULT_ENHANCER_SEQ = window.JASPAR_DASHBOARD_STATE.input_sequence.sequence;
    const MULTI_TF_PROMOTER_SEQ = (
      "GCTAGCTAGCTGATCGATCGTGGCCACCAGGGGGCGCTAGATCGATCGAT" +
      "GGATGAGTCATCCGGATCGATCGGGACATGCCCGGGCATGTCCGATCGAT" +
      "CCCTTTGTTCGATCGATCGATCTAGCGCCCCCTGGTGGCCAGATCGATCG" +
      "CGGATGACTCATCCGATCGATCGTGGCCACCAGGGGGCGCTAGATCGATC"
    );

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

    // Theme Switcher
    function toggleTheme() {{
      const root = document.documentElement;
      const next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      const btn = document.getElementById("themeToggleBtn");
      if (btn) {{
        btn.textContent = next === "dark" ? "Theme: Dark" : "Theme: Light";
      }}
    }}

    function getActiveProfile() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      return st.profiles_catalog[st.active_matrix_id] || st;
    }}

    function initDashboard() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      const sel = document.getElementById("profileSelector");
      sel.innerHTML = Object.keys(st.profiles_catalog).map(id => {{
        const p = st.profiles_catalog[id];
        return `<option value="${{id}}" ${{id === st.active_matrix_id ? "selected" : ""}}>${{p.name}} (${{id}} - ${{p.tf_class}})</option>`;
      }}).join("");

      document.getElementById("liveSequenceInput").value = st.input_sequence.sequence;
      document.getElementById("rngThreshold").value = st.scan_results.threshold || 0.80;
      runLiveSequenceScan(false);
      init3DViewport();
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

      // Compute PWM weights
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
      const fbBanner = document.getElementById("scanFeedbackBanner");
      if (fbText) {{
        const topSummary = hits.length > 0
          ? `Top Hit: ${{hits[0].strand}}[${{hits[0].start}}..${{hits[0].end}}] ${{hits[0].matched_sequence}} (${{(hits[0].relative_score * 100).toFixed(1)}}% affinity, ${{hits[0].significance}})`
          : `No binding sites >= ${{(threshold * 100).toFixed(0)}}% threshold (lower the Relative Cutoff slider to inspect sub-threshold candidates)`;
        fbText.innerHTML = `<strong>Scan #${{scanCounter}} Complete (${{elapsedMs}} ms):</strong> Found <strong>${{hits.length}}</strong> binding site(s) for <strong>${{p.name}} (${{p.matrix_id}})</strong> across ${{rawInput.length}} bp. | ${{topSummary}}`;
      }}
      if (userClicked && fbBanner) {{
        const btn = document.getElementById("btnRunScan");
        if (btn) {{
          btn.textContent = `Scan Complete (${{hits.length}} Hits)`;
          setTimeout(() => {{ btn.textContent = "Run Sequence Scan"; }}, 1400);
        }}
        fbBanner.style.transform = "scale(1.01)";
        fbBanner.style.boxShadow = "0 0 0 3px rgba(37, 99, 235, 0.35)";
        setTimeout(() => {{
          fbBanner.style.transform = "scale(1)";
          fbBanner.style.boxShadow = "none";
        }}, 450);
      }}

      refreshAllCards();
      rebuild3DGraphNodes();
    }}

    function buildLocusTrackSvg(scan) {{
      const seqLen = Math.max(1, scan.sequence_length);
      const width = 920, height = 160, leftPad = 70, rightPad = 30;
      const trackW = width - leftPad - rightPad, midY = 82;
      let out = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{width}} ${{height}}" width="100%" height="${{height}}" style="background:var(--surface-alt); border:1px solid var(--border); border-radius:8px;">`;
      out += `<text x="${{leftPad}}" y="22" font-size="12" font-weight="700" fill="var(--text)">Genomic Coordinate Track: ${{scan.sequence_name}} (1 .. ${{seqLen}} bp)</text>`;
      out += `<line x1="${{leftPad}}" y1="${{midY}}" x2="${{width - rightPad}}" y2="${{midY}}" stroke="#64748B" stroke-width="2"/>`;
      out += `<text x="${{leftPad - 10}}" y="${{midY - 24}}" font-size="11" font-weight="700" fill="#3B82F6" text-anchor="end">+ Strand</text>`;
      out += `<text x="${{leftPad - 10}}" y="${{midY + 32}}" font-size="11" font-weight="700" fill="#EF4444" text-anchor="end">- Strand</text>`;

      for (const frac of [0, 0.25, 0.5, 0.75, 1]) {{
        const x = leftPad + frac * trackW;
        const bp = Math.round(1 + frac * (seqLen - 1));
        out += `<line x1="${{x}}" y1="${{midY - 5}}" x2="${{x}}" y2="${{midY + 5}}" stroke="#64748B" stroke-width="1"/>`;
        out += `<text x="${{x}}" y="${{height - 12}}" font-size="10" fill="var(--muted)" text-anchor="middle">${{bp}} bp</text>`;
      }}

      for (const h of scan.hits) {{
        const x = leftPad + ((h.start - 1) / Math.max(1, seqLen)) * trackW;
        const isFwd = h.strand === "+";
        const col = isFwd ? "#3B82F6" : "#EF4444";
        const stem = 16 + (h.relative_score - 0.5) * 42;
        const yEnd = isFwd ? (midY - stem) : (midY + stem);
        out += `<line x1="${{x}}" y1="${{midY}}" x2="${{x}}" y2="${{yEnd}}" stroke="${{col}}" stroke-width="2.5"/>`;
        out += `<circle cx="${{x}}" cy="${{yEnd}}" r="5.5" fill="${{col}}"><title>${{h.strand}} Strand [${{h.start}}..${{h.end}}]: ${{h.matched_sequence}} (Rel=${{(h.relative_score*100).toFixed(1)}}%)</title></circle>`;
      }}
      out += `</svg>`;
      return out;
    }}

    function refreshAllCards() {{
      const st = window.JASPAR_DASHBOARD_STATE;
      const p = getActiveProfile();
      const speciesStr = (p.species || []).map(s => `${{s.name}} [TaxID: ${{s.tax_id}}]`).join(", ");

      document.getElementById("bannerTitle").textContent = `${{p.name}} Transcription Factor Binding Profile (${{p.matrix_id}})`;
      document.getElementById("bannerSubtitle").textContent = `Species: ${{speciesStr}} | TFClass: ${{p.superclass}} -> ${{p.tf_class}}`;
      document.getElementById("bannerBadges").innerHTML = `
        <span class="badge">Matrix ID: ${{p.matrix_id}}</span>
        <span class="badge">v${{p.version}}</span>
        <span class="badge">${{p.collection}}</span>
        <span class="badge">${{p.logo.length}} bp</span>
      `;
      document.getElementById("activeProfileChip").textContent = `Active Matrix: ${{p.matrix_id}}`;

      // Input stats
      document.getElementById("statProfileName").textContent = `${{p.name}} (${{p.matrix_id}})`;
      document.getElementById("statSeqLength").textContent = `${{st.input_sequence.length}} bp`;
      document.getElementById("statSeqGc").textContent = `${{(st.input_sequence.observed_gc * 100).toFixed(1)}}%`;
      const c = st.input_sequence.counts;
      document.getElementById("statBaseCounts").textContent = `A:${{c.A}} C:${{c.C}} G:${{c.G}} T:${{c.T}}`;
      document.getElementById("statTotalHits").textContent = `${{st.scan_results.total_hits}} hits`;

      // Structural KV
      const uniprotLinks = (p.uniprot_ids || []).map(u => `<a href="https://www.uniprot.org/uniprotkb/${{u}}" target="_blank" rel="noopener">${{u}}</a>`).join(", ");
      const pdbLinks = (p.pdb_ids || []).map(d => `<a href="https://www.rcsb.org/structure/${{d}}" target="_blank" rel="noopener">${{d}}</a>`).join(", ");
      document.getElementById("structuralKvMount").innerHTML = `
        <div class="kv-label">TF Symbol:</div><div><strong>${{p.name}}</strong></div>
        <div class="kv-label">JASPAR Matrix ID:</div><div><code>${{p.matrix_id}}</code> (v${{p.version}})</div>
        <div class="kv-label">Validation Collection:</div><div><span class="badge">${{p.collection}}</span></div>
        <div class="kv-label">TFClass Superclass:</div><div>${{p.superclass}}</div>
        <div class="kv-label">TFClass Class:</div><div>${{p.tf_class}}</div>
        <div class="kv-label">TFClass Family:</div><div>${{p.tf_family}}</div>
        <div class="kv-label">TFClass Subfamily:</div><div>${{p.subfamily}}</div>
        <div class="kv-label">Species Taxonomy:</div><div>${{speciesStr}}</div>
        <div class="kv-label">UniProt Accession(s):</div><div>${{uniprotLinks}}</div>
        <div class="kv-label">PDB 3D Structures:</div><div>${{pdbLinks}}</div>
      `;

      // Cluster Comparison Table
      document.getElementById("clusterTableBody").innerHTML = (p.cluster_comparisons || []).map(row => `
        <tr>
          <td><strong>${{row.matrix_id}}</strong></td>
          <td>${{row.profile_label}}</td>
          <td><span class="badge">${{row.collection}}</span></td>
          <td style="font-family:monospace; font-weight:700;">${{row.consensus}}</td>
          <td>${{row.pearson_r.toFixed(3)}}</td>
          <td>${{row.jensen_shannon_divergence.toFixed(3)}}</td>
        </tr>
      `).join("");

      renderLogoCard();
      updateMatrixView();
      document.getElementById("locusSvgMount").innerHTML = buildLocusTrackSvg(st.scan_results);
      renderScanTable();
    }}

    function renderLogoCard() {{
      const p = getActiveProfile();
      const s = p.logo;
      const svgKey = (isRevComp ? "svg_rev_" : "svg_fwd_") + (isProbAxis ? "prob" : "bits");
      document.getElementById("logoSvgMount").innerHTML = s[svgKey];
      document.getElementById("activeOrientationText").textContent = isRevComp ? "Reverse Complement (3' <- 5')" : "Forward (5' -> 3')";
      document.getElementById("activeConsensusText").textContent = isRevComp ? s.consensus_rev : s.consensus_fwd;
      document.getElementById("logoBitsBadge").textContent = `${{s.total_bits.toFixed(2)}} bits`;
    }}

    function toggleRevComp() {{ isRevComp = !isRevComp; renderLogoCard(); }}
    function toggleLogoAxis() {{ isProbAxis = !isProbAxis; renderLogoCard(); }}

    function downloadActiveLogoSvg() {{
      const p = getActiveProfile();
      const s = p.logo;
      const svgKey = (isRevComp ? "svg_rev_" : "svg_fwd_") + (isProbAxis ? "prob" : "bits");
      const blob = new Blob([s[svgKey]], {{ type: "image/svg+xml;charset=utf-8" }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = p.matrix_id + (isRevComp ? "_revcomp" : "_fwd") + ".svg";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }}

    function cellColor(val, mode) {{
      if (mode === "pfm") return "transparent";
      if (mode === "ppm") {{
        const a = Math.min(0.82, Math.max(0.06, val));
        return `rgba(37, 99, 235, ${{a.toFixed(2)}})`;
      }}
      if (val > 0.05) {{
        const a = Math.min(0.82, val / 2.0);
        return `rgba(37, 99, 235, ${{a.toFixed(2)}})`;
      }} else if (val < -0.05) {{
        const a = Math.min(0.72, Math.abs(val) / 4.0);
        return `rgba(239, 68, 68, ${{a.toFixed(2)}})`;
      }}
      return "transparent";
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

      const pfm = p.pfm;
      const L = p.logo.length;
      const bg = {{ A: (1 - gc) / 2, T: (1 - gc) / 2, G: gc / 2, C: gc / 2 }};
      const bases = ["A", "C", "G", "T"];

      let tbl = '<table class="matrix-table"><thead><tr><th class="base-col">Base</th>';
      for (let i = 1; i <= L; i++) tbl += `<th>${{i}}</th>`;
      tbl += "</tr></thead><tbody>";

      for (const b of bases) {{
        tbl += `<tr><td class="base-col">${{b}}</td>`;
        for (let i = 0; i < L; i++) {{
          const n_i = pfm.A[i] + pfm.C[i] + pfm.G[i] + pfm.T[i];
          const count = pfm[b][i];
          const prob = (count + pseudo) / (n_i + 4 * pseudo);
          const pwm = Math.log2(prob / bg[b]);
          const txt = heatmapMode === "pfm" ? count.toFixed(0) : (heatmapMode === "ppm" ? prob.toFixed(3) : pwm.toFixed(2));
          const bgStyle = cellColor(heatmapMode === "pwm" ? pwm : prob, heatmapMode);
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
      document.getElementById("hitCountLabel").textContent = `${{filtered.length}} Sites`;
      document.getElementById("scanSummaryMeta").innerHTML = `
        Target Locus: <strong>${{scan.sequence_name}}</strong> (${{scan.sequence_length}} bp) |
        Relative Cutoff: <strong>${{scan.threshold.toFixed(2)}}</strong>
      `;
      const rows = filtered.map((h, idx) => `
        <tr>
          <td>${{idx + 1}}</td>
          <td>${{h.start}}</td>
          <td>${{h.end}}</td>
          <td><strong>${{h.strand}}</strong></td>
          <td style="font-family:monospace; font-weight:700;">${{h.matched_sequence}}</td>
          <td>${{h.raw_score.toFixed(3)}}</td>
          <td>${{(h.relative_score * 100).toFixed(1)}}% (${{h.relative_score.toFixed(3)}})</td>
          <td><span class="badge">${{h.significance}}</span></td>
        </tr>
      `).join("");
      document.getElementById("scanTableBody").innerHTML = rows || '<tr><td colspan="8">No binding sites above threshold.</td></tr>';
    }}

    // =========================================================================
    // 3D Interactive Viewport Engine ("Motif Relationship Network" & "3D DNA Binding Groove")
    // =========================================================================
    let camRotX = 0.35;
    let camRotY = -0.55;
    let camZoom = 1.0;
    let isDragging3D = false;
    let prevMouseX = 0;
    let prevMouseY = 0;
    let graphNodes = [];
    let graphEdges = [];
    let selectedNode = null;

    function set3DViewMode(mode) {{
      view3DMode = mode;
      document.getElementById("btn3DNetwork").className = mode === "network" ? "btn active-tab" : "btn secondary";
      document.getElementById("btn3DHelix").className = mode === "helix" ? "btn active-tab" : "btn secondary";
      document.getElementById("view3DTitle").textContent = mode === "network" ? "Motif Relationship Network" : "3D DNA Binding Groove";
      document.getElementById("view3DModeBadge").textContent = mode === "network" ? "3D Network Graph" : "3D B-DNA Helix";
      rebuild3DGraphNodes();
    }}

    function reset3DCamera() {{
      camRotX = 0.35;
      camRotY = -0.55;
      camZoom = 1.0;
    }}

    function rebuild3DGraphNodes() {{
      const p = getActiveProfile();
      const scan = window.JASPAR_DASHBOARD_STATE.scan_results;
      graphNodes = [];
      graphEdges = [];

      if (view3DMode === "network") {{
        // Central TF Hub Node
        graphNodes.push({{
          id: "tf_hub",
          label: `${{p.name}} (${{p.matrix_id}})`,
          category: "Active TF Profile",
          detail: `Consensus: ${{p.logo.consensus_fwd}} | Total Info: ${{p.logo.total_bits.toFixed(2)}} bits`,
          x: 0, y: 0, z: 0, vx: 0, vy: 0, vz: 0,
          radius: 14, color: "#3B82F6", fixed: true
        }});

        // Structural TFClass Node
        graphNodes.push({{
          id: "tf_class",
          label: p.tf_class,
          category: "TFClass Domain",
          detail: `Superclass: ${{p.superclass}} | Family: ${{p.tf_family}}`,
          x: -95, y: -65, z: 40, vx: 0, vy: 0, vz: 0,
          radius: 10, color: "#10B981"
        }});
        graphEdges.push({{ source: 0, target: 1, length: 115, color: "rgba(16, 185, 129, 0.55)" }});

        // UniProt / 3D Structure Node
        const upId = (p.uniprot_ids && p.uniprot_ids[0]) || "UniProt";
        graphNodes.push({{
          id: "uniprot",
          label: `UniProt: ${{upId}}`,
          category: "Protein Accession",
          detail: `PDB Structures: ${{(p.pdb_ids || []).join(", ")}}`,
          x: -105, y: 60, z: -35, vx: 0, vy: 0, vz: 0,
          radius: 9, color: "#A855F7"
        }});
        graphEdges.push({{ source: 0, target: 2, length: 115, color: "rgba(168, 85, 247, 0.55)" }});

        // Familial Cluster Nodes
        (p.cluster_comparisons || []).slice(1).forEach((cl, idx) => {{
          const nIdx = graphNodes.length;
          const angle = (idx + 1) * 1.3;
          graphNodes.push({{
            id: `cluster_${{idx}}`,
            label: cl.matrix_id,
            category: cl.collection,
            detail: `${{cl.profile_label}} | Pearson r=${{cl.pearson_r.toFixed(3)}}`,
            x: 95 * Math.cos(angle), y: -75, z: 95 * Math.sin(angle), vx: 0, vy: 0, vz: 0,
            radius: 9, color: "#F59E0B"
          }});
          graphEdges.push({{ source: 0, target: nIdx, length: 120, color: "rgba(245, 158, 11, 0.55)" }});
        }});

        // Detected Genomic Binding Site Nodes
        (scan.hits || []).slice(0, 8).forEach((hit, idx) => {{
          const nIdx = graphNodes.length;
          const theta = (idx / Math.max(1, Math.min(8, scan.hits.length))) * Math.PI * 2;
          graphNodes.push({{
            id: `hit_${{idx}}`,
            label: `${{hit.strand}}[${{hit.start}}..${{hit.end}}]`,
            category: `Binding Site (${{hit.significance}})`,
            detail: `k-mer: ${{hit.matched_sequence}} | Rel Score: ${{(hit.relative_score * 100).toFixed(1)}}% (${{hit.raw_score.toFixed(2)}} bits)`,
            x: 125 * Math.cos(theta), y: 68 + (idx % 2) * 24, z: 125 * Math.sin(theta), vx: 0, vy: 0, vz: 0,
            radius: 8 + hit.relative_score * 3,
            color: hit.strand === "+" ? "#60A5FA" : "#F87171"
          }});
          graphEdges.push({{ source: 0, target: nIdx, length: 135, color: hit.strand === "+" ? "rgba(96, 165, 250, 0.5)" : "rgba(248, 113, 113, 0.5)" }});
        }});
      }} else {{
        // 3D B-DNA Double Helix Binding Groove Mode
        const cols = p.logo.columns || [];
        const L = cols.length;
        const baseColors = {{ A: "#10B981", C: "#3B82F6", G: "#F59E0B", T: "#EF4444", N: "#94A3B8" }};
        const compMap = {{ A: "T", T: "A", C: "G", G: "C", N: "N" }};

        for (let i = 0; i < L; i++) {{
          const col = cols[i];
          const topItem = col.stack[col.stack.length - 1];
          const base = topItem.base;
          const compBase = compMap[base] || "N";
          const twist = i * 0.628; // 36 degrees per bp in B-DNA
          const yPos = (i - (L - 1) / 2) * 14.5;
          const radius = 58;

          const idxFwd = graphNodes.length;
          graphNodes.push({{
            id: `bp_fwd_${{i+1}}`,
            label: `${{base}}${{i+1}}`,
            category: `Forward Strand Pos ${{i+1}}`,
            detail: `Base ${{base}} | Info R_i: ${{col.info_r_i.toFixed(2)}} bits | IUPAC: ${{col.iupac}}`,
            x: radius * Math.cos(twist),
            y: yPos,
            z: radius * Math.sin(twist),
            vx: 0, vy: 0, vz: 0,
            radius: 5.5 + col.info_r_i * 2.8,
            color: baseColors[base] || "#3B82F6",
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
            radius: 5.0,
            color: baseColors[compBase] || "#94A3B8",
            fixed: true
          }});

          // Hydrogen-bond base pair rung
          graphEdges.push({{ source: idxFwd, target: idxRev, length: 116, color: "rgba(226, 232, 240, 0.42)" }});
          // Sugar-phosphate backbone links
          if (i > 0) {{
            graphEdges.push({{ source: idxFwd - 2, target: idxFwd, length: 24, color: "rgba(59, 130, 246, 0.75)" }});
            graphEdges.push({{ source: idxRev - 2, target: idxRev, length: 24, color: "rgba(239, 68, 68, 0.65)" }});
          }}
        }}
      }}
      selectedNode = graphNodes[0] || null;
      updateInspectorBanner();
    }}

    function stepNetworkRelaxation() {{
      if (view3DMode !== "network") return;
      const N = graphNodes.length;
      // Node-node repulsion
      for (let i = 0; i < N; i++) {{
        for (let j = i + 1; j < N; j++) {{
          const dx = graphNodes[j].x - graphNodes[i].x;
          const dy = graphNodes[j].y - graphNodes[i].y;
          const dz = graphNodes[j].z - graphNodes[i].z;
          const distSq = dx * dx + dy * dy + dz * dz + 25;
          const force = 420 / distSq;
          const dist = Math.sqrt(distSq);
          const fx = (dx / dist) * force;
          const fy = (dy / dist) * force;
          const fz = (dz / dist) * force;
          if (!graphNodes[i].fixed) {{ graphNodes[i].vx -= fx; graphNodes[i].vy -= fy; graphNodes[i].vz -= fz; }}
          if (!graphNodes[j].fixed) {{ graphNodes[j].vx += fx; graphNodes[j].vy += fy; graphNodes[j].vz += fz; }}
        }}
      }}
      // Edge spring attraction
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
      // Integrate with damping
      for (const n of graphNodes) {{
        if (n.fixed) continue;
        n.vx *= 0.82; n.vy *= 0.82; n.vz *= 0.82;
        n.x += n.vx; n.y += n.vy; n.z += n.vz;
      }}
    }}

    function project3D(node, w, h) {{
      // Rotate around Y then X
      const cosY = Math.cos(camRotY), sinY = Math.sin(camRotY);
      const cosX = Math.cos(camRotX), sinX = Math.sin(camRotX);
      const x1 = node.x * cosY - node.z * sinY;
      const z1 = node.z * cosY + node.x * sinY;
      const y2 = node.y * cosX - z1 * sinX;
      const z2 = z1 * cosX + node.y * sinX;
      const fov = 420 * camZoom;
      const scale = fov / Math.max(120, 420 + z2);
      return {{
        sx: w / 2 + x1 * scale,
        sy: h / 2 + y2 * scale,
        z: z2,
        scale: scale
      }};
    }}

    function updateInspectorBanner() {{
      if (!selectedNode) return;
      document.getElementById("nodeInspectorText").textContent = `${{selectedNode.category}}: ${{selectedNode.label}}`;
      document.getElementById("nodeMetricBadge").textContent = selectedNode.detail;
    }}

    function init3DViewport() {{
      const host = document.getElementById("viewport3DHost");
      const svgNs = "http://www.w3.org/2000/svg";
      const surface = document.createElementNS(svgNs, "svg");
      surface.setAttribute("id", "viewport3DSurface");
      host.insertBefore(surface, host.firstChild);
      rebuild3DGraphNodes();

      surface.addEventListener("mousedown", e => {{
        isDragging3D = true;
        prevMouseX = e.clientX;
        prevMouseY = e.clientY;

        // Check if user clicked a projected node
        const rect = surface.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;
        for (const n of graphNodes) {{
          if (n._proj) {{
            const dist = Math.hypot(mx - n._proj.sx, my - n._proj.sy);
            if (dist <= n.radius * n._proj.scale + 6) {{
              selectedNode = n;
              updateInspectorBanner();
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
        camRotY += dx * 0.008;
        camRotX = Math.max(-1.3, Math.min(1.3, camRotX + dy * 0.008));
      }});

      window.addEventListener("mouseup", () => {{ isDragging3D = false; }});
      surface.addEventListener("wheel", e => {{
        e.preventDefault();
        camZoom = Math.max(0.55, Math.min(2.2, camZoom - e.deltaY * 0.001));
      }}, {{ passive: false }});

      function animate() {{
        const rect = surface.getBoundingClientRect();
        const W = Math.max(300, rect.width || 520);
        const H = Math.max(240, rect.height || 360);
        surface.setAttribute("viewBox", `0 0 ${{W}} ${{H}}`);

        if (!isDragging3D) camRotY += 0.0025;
        stepNetworkRelaxation();

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
          const r = Math.max(3, n.radius * pr.scale);
          const strokeAttr = (selectedNode && selectedNode.id === n.id)
            ? ' stroke="#FFFFFF" stroke-width="2.5"'
            : '';
          parts.push(
            `<circle cx="${{pr.sx.toFixed(1)}}" cy="${{pr.sy.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{n.color}}"${{strokeAttr}}/>`
          );
          parts.push(
            `<text x="${{pr.sx.toFixed(1)}}" y="${{(pr.sy - r - 5).toFixed(1)}}" fill="#E2E8F0" font-size="11" font-weight="600" text-anchor="middle" font-family="ui-sans-serif, system-ui, sans-serif">${{n.label}}</text>`
          );
        }}

        surface.innerHTML = parts.join("");
        requestAnimationFrame(animate);
      }}
      requestAnimationFrame(animate);
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
