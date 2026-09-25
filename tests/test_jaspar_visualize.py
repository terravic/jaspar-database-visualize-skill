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
"""Unit, Smoke, Mathematical, and End-to-End System Tests for jaspar-database-visualize."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Dict, List

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True

import pytest

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = WORKSPACE_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from generate_dashboard import compile_dashboard_html
from jaspar_api import (
    compute_ppm_and_pwm,
    format_motif_output,
    load_mock_fixture,
    verify_workspace_license,
)
from visualize_logo import (
    compute_logo_metrics,
    derive_iupac_symbol,
    render_logo_svg,
    reverse_complement_pfm,
)
from visualize_matrix import build_matrix_analysis
from visualize_scan import scan_dna_sequence


def test_legal_compliance_and_apache_headers() -> None:
    """Verifies workspace license guard, SKILL_LICENSES.md, and Apache 2.0 headers."""
    lic_path = verify_workspace_license(quiet=True)
    assert lic_path.exists()
    lic_text = lic_path.read_text(encoding="utf-8")
    assert "Creative Commons Attribution 4.0 International (CC BY 4.0)" in lic_text
    assert "Timestamp:" in lic_text

    skill_licenses = WORKSPACE_ROOT / "SKILL_LICENSES.md"
    assert skill_licenses.exists()
    assert "jaspar-database-visualize" in skill_licenses.read_text(encoding="utf-8")

    expected_scripts = [
        "jaspar_api.py",
        "visualize_logo.py",
        "visualize_matrix.py",
        "visualize_scan.py",
        "generate_dashboard.py",
        "run_pipeline.py",
    ]
    for script_name in expected_scripts:
        script_path = SCRIPTS_DIR / script_name
        assert script_path.exists(), f"Missing script {script_path}"
        content = script_path.read_text(encoding="utf-8")
        assert "Licensed under the Apache License, Version 2.0" in content
        assert "Copyright 2026 Project Contributors" in content


def test_ppm_and_pwm_mathematics() -> None:
    """Verifies PPM pseudocount normalization and PWM log2-odds scoring."""
    pfm: Dict[str, List[float]] = {
        "A": [100, 0],
        "C": [0, 100],
        "G": [0, 0],
        "T": [0, 0],
    }
    res = compute_ppm_and_pwm(pfm, pseudocount=0.8, gc_content=0.50)
    assert res["length"] == 2
    # Column 1: N_1 = 100, denom = 103.2, PPM[A, 0] = 100.8 / 103.2
    expected_p_a = 100.8 / 103.2
    assert math.isclose(res["ppm"]["A"][0], expected_p_a, rel_tol=1e-5)
    # Uniform background = 0.25 -> PWM[A, 0] = log2(expected_p_a / 0.25)
    expected_w_a = math.log2(expected_p_a / 0.25)
    assert math.isclose(res["pwm"]["A"][0], expected_w_a, rel_tol=1e-5)

    with pytest.raises(ValueError):
        compute_ppm_and_pwm(pfm, pseudocount=0.8, gc_content=0.0)


def test_shannon_entropy_and_iupac_consensus() -> None:
    """Verifies Shannon entropy H_i, Schneider correction e_n, R_i, and IUPAC symbols."""
    assert derive_iupac_symbol({"A": 0.80, "C": 0.10, "G": 0.05, "T": 0.05}) == "A"
    assert derive_iupac_symbol({"A": 0.45, "G": 0.40, "C": 0.10, "T": 0.05}) == "R"
    assert derive_iupac_symbol({"C": 0.42, "T": 0.40, "A": 0.10, "G": 0.08}) == "Y"
    assert derive_iupac_symbol({"A": 0.30, "G": 0.30, "T": 0.30, "C": 0.10}) == "D"
    assert derive_iupac_symbol({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}) == "N"

    ctcf_data = load_mock_fixture("MA0139.1")
    fwd = compute_logo_metrics(ctcf_data["pfm"], pseudocount=0.8, reverse_complement=False)
    rev = compute_logo_metrics(ctcf_data["pfm"], pseudocount=0.8, reverse_complement=True)

    assert fwd["length"] == 19
    assert rev["length"] == 19
    assert math.isclose(fwd["total_bits"], rev["total_bits"], rel_tol=1e-5)
    for col in fwd["columns"]:
        assert 0.0 <= col["info_r_i"] <= 2.0
        # Verify stack is sorted ascending by bits
        bits_list = [item["bits"] for item in col["stack"]]
        assert bits_list == sorted(bits_list)

    svg_doc = render_logo_svg(fwd, title="CTCF Test Logo", mode="bits")
    assert "<svg" in svg_doc
    assert "<path" in svg_doc
    assert "#10B981" in svg_doc  # A green
    assert "#3B82F6" in svg_doc  # C blue


def test_dual_strand_sequence_scanner_and_edge_cases() -> None:
    """Verifies forward (+) and reverse (-) strand binding site detection and boundary guards."""
    ctcf_data = load_mock_fixture("MA0139.1")
    fasta_path = WORKSPACE_ROOT / "sample_data" / "sample_enhancer.fasta"
    lines = [ln.strip() for ln in fasta_path.read_text(encoding="utf-8").splitlines() if not ln.startswith(">")]
    enhancer_seq = "".join(lines)
    assert len(enhancer_seq) == 1000

    scan_res = scan_dna_sequence(
        matrix_data=ctcf_data,
        sequence=enhancer_seq,
        sequence_name="SYNTHETIC_HUMAN_ENHANCER_1000BP",
        threshold=0.80,
    )
    assert scan_res["total_hits"] >= 3
    strands = {h["strand"] for h in scan_res["hits"]}
    assert "+" in strands
    assert "-" in strands

    # Verify short sequence handling
    short_res = scan_dna_sequence(ctcf_data, "ACGT", threshold=0.80)
    assert short_res["total_hits"] == 0
    assert len(short_res["warnings"]) == 1

    # Verify >100kb sequence guardrail
    with pytest.raises(ValueError):
        scan_dna_sequence(ctcf_data, "A" * 100005, threshold=0.80)


def test_dashboard_generation_and_mandatory_citations() -> None:
    """Verifies standalone HTML5 dashboard contains all required cards, citations, and disclaimers."""
    ctcf_data = load_mock_fixture("MA0139.1")
    html_doc = compile_dashboard_html(motif_data=ctcf_data)

    assert "<!DOCTYPE html>" in html_doc
    assert "window.JASPAR_DASHBOARD_STATE" in html_doc
    assert "Ovek Baydar, D. et al." in html_doc
    assert "Rauluseviciute, I. et al." in html_doc
    assert "doi:10.1093/nar/gkaf1209" in html_doc
    assert "Creative Commons Attribution 4.0 International (CC BY 4.0)" in html_doc
    assert "Predicted binding sites must be experimentally validated via ChIP-seq, CUT&RUN, or reporter assays." in html_doc


def test_cli_acceptance_suite() -> None:
    """Executes all 6 acceptance CLI workflows from Section 9.1 in a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        resolve_out = tmp / "ctcf_resolve.json"
        motif_out = tmp / "ctcf_motif.json"
        logo_svg = tmp / "ctcf_logo.svg"
        logo_html = tmp / "ctcf_logo.html"
        matrix_html = tmp / "ctcf_matrix.html"
        hits_json = tmp / "ctcf_hits.json"
        dash_html = tmp / "ctcf_dashboard.html"
        pipe_dir = tmp / "ctcf_pipeline_output"
        fasta_file = WORKSPACE_ROOT / "sample_data" / "sample_enhancer.fasta"

        cmds = [
            [
                sys.executable,
                str(SCRIPTS_DIR / "jaspar_api.py"),
                "resolve_tf_id",
                "--name",
                "CTCF",
                "--tax-id",
                "9606",
                "--output",
                str(resolve_out),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "jaspar_api.py"),
                "get_tf_motif",
                "--matrix-id",
                "MA0139.1",
                "--output",
                str(motif_out),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "visualize_logo.py"),
                "--matrix-json",
                str(motif_out),
                "--output-svg",
                str(logo_svg),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "visualize_logo.py"),
                "--matrix-json",
                str(motif_out),
                "--output-html",
                str(logo_html),
                "--reverse-complement",
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "visualize_matrix.py"),
                "--matrix-json",
                str(motif_out),
                "--output-html",
                str(matrix_html),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "visualize_scan.py"),
                "--matrix-json",
                str(motif_out),
                "--fasta",
                str(fasta_file),
                "--threshold",
                "0.80",
                "--output-json",
                str(hits_json),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "generate_dashboard.py"),
                "--matrix-id",
                "MA0139.1",
                "--motif-json",
                str(motif_out),
                "--scan-json",
                str(hits_json),
                "--output",
                str(dash_html),
                "--mock",
            ],
            [
                sys.executable,
                str(SCRIPTS_DIR / "run_pipeline.py"),
                "--name",
                "CTCF",
                "--tax-id",
                "9606",
                "--fasta",
                str(fasta_file),
                "--output-dir",
                str(pipe_dir),
                "--mock",
            ],
        ]

        for cmd in cmds:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            assert proc.returncode == 0, f"Command failed: {' '.join(cmd)}\nSTDERR: {proc.stderr}"

        assert resolve_out.exists()
        assert motif_out.exists()
        assert logo_svg.exists()
        assert logo_html.exists()
        assert matrix_html.exists()
        assert hits_json.exists()
        assert dash_html.exists()
        assert (pipe_dir / "dashboard.html").exists()


def test_custom_user_data_pipeline_and_dashboard_reflection() -> None:
    """Verifies that custom user DNA sequences and TFs outside sample_data/ are processed and reflected in dashboard.html."""
    custom_user_seq = "GCTAGCTAGCGACCACGTGGTCGATCGATCGGACCACGTGGTCGATCGATCG"
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir) / "myc_custom_output"
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "run_pipeline.py"),
            "--name",
            "MYC",
            "--tax-id",
            "9606",
            "--sequence",
            custom_user_seq,
            "--threshold",
            "0.80",
            "--output-dir",
            str(out_dir),
            "--mock",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        assert proc.returncode == 0, f"Pipeline failed on custom user data: {proc.stderr}"

        dash_file = out_dir / "dashboard.html"
        assert dash_file.exists()
        dash_html = dash_file.read_text(encoding="utf-8")
        # Verify the processed MYC matrix (MA0147.3) and custom user sequence are embedded in dashboard.html
        assert "MA0147.3" in dash_html
        assert "MYC" in dash_html
        assert custom_user_seq in dash_html

