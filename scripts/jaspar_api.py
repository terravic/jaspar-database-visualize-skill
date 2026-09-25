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
"""JASPAR REST API v1 Client and Deterministic Offline Mock Engine.

Supports all baseline jaspar-database subcommands (resolve_tf_id, get_tf_motif,
get_tf_metadata, get_tf_pwm, infer_from_sequence, get_tffm) while introducing:
1. Automated workspace license verification (.licenses/jaspar_database_visualize_LICENSE.txt)
2. Persistent file output (--output <file>) across all subcommands
3. Deterministic offline mock mode (--mock) backed by local synthetic fixtures
4. Multi-format matrix serialization (json, jsonp, jaspar, meme, transfac, pfm, yaml)
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional
import urllib.error
import urllib.parse
import urllib.request

JASPAR_API_BASE = "https://jaspar.elixir.no/api/v1"
MAX_STDOUT_CHARS = 50000

ACADEMIC_CITATION_NOTICE = (
    "Citations:\n"
    "  1. Ovek Baydar, D. et al. JASPAR 2026: expansion of transcription factor binding "
    "profiles and integration of deep learning models. Nucleic Acids Research 54, D184-D193 (2026). "
    "doi:10.1093/nar/gkaf1209.\n"
    "  2. Rauluseviciute, I. et al. JASPAR 2024: 20th anniversary of the open-access database "
    "of transcription factor binding profiles. Nucleic Acids Research 52, D174-D182 (2024). "
    "doi:10.1093/nar/gkad1059.\n"
    "License: Creative Commons Attribution 4.0 International (CC BY 4.0)."
)

SCIENTIFIC_DISCLAIMER = (
    "JASPAR position frequency matrices represent in vitro or statistical models of sequence "
    "affinity. In vivo transcription factor binding depends strongly on chromatin accessibility, "
    "nucleosome positioning, DNA methylation, and cofactor interactions. Predicted binding sites "
    "must be experimentally validated via ChIP-seq, CUT&RUN, or reporter assays."
)

MOCK_MATRIX_MAP: Dict[str, str] = {
    "MA0139.1": "ctcf_ma0139_1.json",
    "MA0139": "ctcf_ma0139_1.json",
    "MA0488.2": "jun_ma0488_2.json",
    "MA0488": "jun_ma0488_2.json",
    "MA0106.3": "tp53_ma0106_3.json",
    "MA0106": "tp53_ma0106_3.json",
    "MA0143.4": "sox2_ma0143_4.json",
    "MA0143": "sox2_ma0143_4.json",
}

MOCK_NAME_MAP: Dict[str, str] = {
    "CTCF": "MA0139.1",
    "JUN": "MA0488.2",
    "AP-1": "MA0488.2",
    "AP1": "MA0488.2",
    "TP53": "MA0106.3",
    "P53": "MA0106.3",
    "SOX2": "MA0143.4",
}


def find_skill_root() -> Path:
    """Locates the project root directory containing sample_data/ and scripts/."""
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent


def find_workspace_root() -> Path:
    """Locates the top-level workspace root containing SKILL_LICENSES.md or .licenses."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "SKILL_LICENSES.md").exists() or (parent / ".licenses").is_dir():
            return parent
    return Path(__file__).resolve().parent.parent


def verify_workspace_license(quiet: bool = False) -> Path:
    """Verifies or creates the workspace license verification file per Section 3.2."""
    workspace_root = find_workspace_root()
    license_dir = workspace_root / ".licenses"
    primary_file = license_dir / "jaspar_database_visualize_LICENSE.txt"
    fallback_file = license_dir / "jaspar_database_LICENSE.txt"

    if primary_file.exists():
        return primary_file
    if fallback_file.exists():
        return fallback_file

    license_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    notice_text = (
        "Notification: Please check the JASPAR terms of use and licensing at "
        "https://jaspar.elixir.no/ and https://jaspar.elixir.no/api/. "
        "JASPAR data is available under Creative Commons Attribution 4.0 International (CC BY 4.0).\n"
        f"Timestamp: {timestamp}\n"
    )
    primary_file.write_text(notice_text, encoding="utf-8")
    if not quiet:
        sys.stderr.write(f"[LICENSE GUARD] Created {primary_file}\n{notice_text}")
    return primary_file


class PoliteHttpClient:
    """Rate-limited HTTP client enforcing 10 QPS and exponential backoff retry."""

    def __init__(self, qps: float = 10.0, max_retries: int = 3) -> None:
        self.min_interval = 1.0 / max(qps, 0.1)
        self.max_retries = max_retries
        self._last_request_ts = 0.0
        self._external_client: Optional[Any] = None
        try:
            from polite_http.http_client import HttpClient  # type: ignore
            self._external_client = HttpClient(qps=qps)
        except Exception:
            self._external_client = None

    def get_json(self, url: str) -> Dict[str, Any]:
        """Fetches JSON payload from target URL with rate limiting and backoff."""
        elapsed = time.monotonic() - self._last_request_ts
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)

        for attempt in range(self.max_retries):
            try:
                self._last_request_ts = time.monotonic()
                req = urllib.request.Request(
                    url,
                    headers={
                        "Accept": "application/json",
                        "User-Agent": "jaspar-database-visualize/1.0 (Apache-2.0; +https://jaspar.elixir.no)",
                    },
                )
                with urllib.request.urlopen(req, timeout=15.0) as resp:
                    raw = resp.read().decode("utf-8")
                    return json.loads(raw)
            except (urllib.error.HTTPError, urllib.error.URLError) as exc:
                if attempt + 1 >= self.max_retries:
                    raise RuntimeError(
                        f"JASPAR API request failed after {self.max_retries} attempts for URL '{url}': {exc}. "
                        "Pass --mock for deterministic offline execution."
                    ) from exc
                time.sleep(0.25 * (2 ** attempt))
        raise RuntimeError(f"Unreachable JASPAR endpoint: {url}")


EXTENDED_CANONICAL_PROFILES: Dict[str, Dict[str, Any]] = {
    "MA0147.3": {
        "matrix_id": "MA0147.3", "base_id": "MA0147", "version": "3", "name": "MYC",
        "collection": "CORE", "tax_group": "vertebrates",
        "species": [{"tax_id": 9606, "name": "Homo sapiens"}],
        "superclass": "Basic domains", "class": ["Basic helix-loop-helix factors (bHLH)"],
        "family": ["bHLH-ZIP factors"], "subfamily": "MYC/MAX heterodimer E-box regulators",
        "type": "ChIP-seq", "uniprot_ids": ["P01106"], "pdb_ids": ["1NKP", "6G6K"],
        "pfm": {
            "A": [120, 85, 40, 955, 15, 20, 945, 25, 30, 115, 140, 160],
            "C": [240, 790, 920, 15, 950, 25, 15, 18, 22, 780, 230, 210],
            "G": [480, 65, 25, 18, 20, 940, 22, 935, 925, 45, 460, 440],
            "T": [160, 60, 15, 12, 15, 15, 18, 22, 23, 60, 170, 190],
        },
    },
    "MA0148.4": {
        "matrix_id": "MA0148.4", "base_id": "MA0148", "version": "4", "name": "FOXA1",
        "collection": "CORE", "tax_group": "vertebrates",
        "species": [{"tax_id": 9606, "name": "Homo sapiens"}],
        "superclass": "Helix-turn-helix domains", "class": ["Fork head / winged helix factors"],
        "family": ["Forkhead box (FOX) factors"], "subfamily": "Pioneer chromatin opening factors",
        "type": "ChIP-seq", "uniprot_ids": ["P55317"], "pdb_ids": ["1VTN", "7VOX"],
        "pfm": {
            "A": [220, 140, 35, 920, 940, 45, 950, 890, 930, 910, 180],
            "C": [210, 160, 40, 25, 20, 30, 15, 35, 25, 30, 640],
            "G": [290, 520, 30, 30, 20, 905, 15, 40, 20, 25, 80],
            "T": [280, 180, 895, 25, 20, 20, 20, 35, 25, 35, 100],
        },
    },
    "MA0035.4": {
        "matrix_id": "MA0035.4", "base_id": "MA0035", "version": "4", "name": "GATA1",
        "collection": "CORE", "tax_group": "vertebrates",
        "species": [{"tax_id": 9606, "name": "Homo sapiens"}],
        "superclass": "Zinc-coordinating DNA-binding domains", "class": ["Other C4 zinc finger-type factors"],
        "family": ["GATA-type zinc fingers"], "subfamily": "Erythroid lineage GATA regulators",
        "type": "ChIP-seq", "uniprot_ids": ["P15976"], "pdb_ids": ["1GNF", "3DFV"],
        "pfm": {
            "A": [210, 180, 40, 950, 20, 940, 935, 85, 160, 220, 190],
            "C": [520, 580, 35, 15, 20, 20, 25, 45, 680, 260, 250],
            "G": [140, 120, 905, 20, 25, 20, 25, 845, 90, 280, 310],
            "T": [130, 120, 20, 15, 935, 20, 15, 25, 70, 240, 250],
        },
    },
}

EXTENDED_NAME_LOOKUP: Dict[str, str] = {
    "MYC": "MA0147.3",
    "MA0147": "MA0147.3",
    "MA0147.3": "MA0147.3",
    "FOXA1": "MA0148.4",
    "MA0148": "MA0148.4",
    "MA0148.4": "MA0148.4",
    "GATA1": "MA0035.4",
    "MA0035": "MA0035.4",
    "MA0035.4": "MA0035.4",
}


def load_mock_fixture(matrix_id_or_name: str) -> Dict[str, Any]:
    """Loads a synthetic matrix JSON fixture from sample_data/ or built-in canonical profiles."""
    skill_root = find_skill_root()
    sample_dir = skill_root / "sample_data"
    key = matrix_id_or_name.strip().upper()
    if key in EXTENDED_NAME_LOOKUP:
        ext_id = EXTENDED_NAME_LOOKUP[key]
        return dict(EXTENDED_CANONICAL_PROFILES[ext_id])
    if key in MOCK_NAME_MAP:
        key = MOCK_NAME_MAP[key]

    filename = MOCK_MATRIX_MAP.get(key)
    if not filename:
        # Check prefix match (e.g., MA0139.2 -> ctcf_ma0139_1.json)
        base_id = key.split(".")[0]
        filename = MOCK_MATRIX_MAP.get(base_id)

    if not filename:
        supported = ", ".join(sorted(set(MOCK_NAME_MAP.values()) | set(EXTENDED_CANONICAL_PROFILES.keys())))
        raise KeyError(
            f"Unknown mock matrix or TF identifier '{matrix_id_or_name}'. "
            f"Supported mock fixtures: {supported} (CTCF, JUN, TP53, SOX2, MYC, FOXA1, GATA1)."
        )

    fixture_path = sample_dir / filename
    if not fixture_path.exists():
        raise FileNotFoundError(f"Mock fixture file missing at {fixture_path}")
    return json.loads(fixture_path.read_text(encoding="utf-8"))


def validate_pfm(pfm: Dict[str, List[float]]) -> int:
    """Validates that a PFM dictionary contains A, C, G, T arrays of equal positive length."""
    required_bases = ("A", "C", "G", "T")
    for base in required_bases:
        if base not in pfm or not isinstance(pfm[base], list):
            raise ValueError(f"PFM dictionary is missing required nucleotide array '{base}'.")
    lengths = {len(pfm[b]) for b in required_bases}
    if len(lengths) != 1:
        raise ValueError(f"PFM nucleotide rows have unequal lengths: { {b: len(pfm[b]) for b in required_bases} }.")
    length = lengths.pop()
    if length < 1:
        raise ValueError("PFM matrix must have at least 1 column position.")
    for base in required_bases:
        for idx, val in enumerate(pfm[base]):
            if float(val) < 0.0:
                raise ValueError(f"Negative count {val} at base {base}, position {idx + 1}.")
    return length


def compute_ppm_and_pwm(
    pfm: Dict[str, List[float]],
    pseudocount: float = 0.8,
    gc_content: float = 0.50,
) -> Dict[str, Any]:
    """Computes normalized Position Probability Matrix (PPM) and log2 Position Weight Matrix (PWM)."""
    if pseudocount <= 0.0 or pseudocount > 100.0:
        raise ValueError(f"Pseudocount must be in (0.0, 100.0], received {pseudocount}.")
    if not (0.01 <= gc_content <= 0.99):
        raise ValueError(f"GC content must be strictly within [0.01, 0.99], received {gc_content}.")

    length = validate_pfm(pfm)
    bg = {
        "A": (1.0 - gc_content) / 2.0,
        "T": (1.0 - gc_content) / 2.0,
        "G": gc_content / 2.0,
        "C": gc_content / 2.0,
    }
    ppm: Dict[str, List[float]] = {b: [] for b in ("A", "C", "G", "T")}
    pwm: Dict[str, List[float]] = {b: [] for b in ("A", "C", "G", "T")}
    col_totals: List[float] = []

    for i in range(length):
        n_i = sum(float(pfm[b][i]) for b in ("A", "C", "G", "T"))
        col_totals.append(n_i)
        denom = n_i + 4.0 * pseudocount
        for b in ("A", "C", "G", "T"):
            prob = (float(pfm[b][i]) + pseudocount) / denom
            weight = math.log2(prob / bg[b])
            ppm[b].append(round(prob, 6))
            pwm[b].append(round(weight, 6))

    return {
        "length": length,
        "pseudocount": pseudocount,
        "gc_content": gc_content,
        "background": bg,
        "column_totals": col_totals,
        "ppm": ppm,
        "pwm": pwm,
    }


def format_motif_output(data: Dict[str, Any], fmt: str) -> str:
    """Converts matrix profile dictionary into requested standard motif format."""
    fmt_clean = fmt.strip().lower()
    matrix_id = data.get("matrix_id", "MA0000.1")
    name = data.get("name", "UNKNOWN")
    pfm = data.get("pfm", {})
    length = validate_pfm(pfm)

    if fmt_clean == "json":
        return json.dumps(data, indent=2)
    if fmt_clean == "jsonp":
        return f"jasparCallback({json.dumps(data, indent=2)});"
    if fmt_clean == "pfm":
        lines = []
        for b in ("A", "C", "G", "T"):
            row_str = " ".join(f"{int(round(v))}" for v in pfm[b])
            lines.append(row_str)
        return "\n".join(lines) + "\n"
    if fmt_clean == "jaspar":
        lines = [f">{matrix_id} {name}"]
        for b in ("A", "C", "G", "T"):
            row_str = " ".join(f"{int(round(v)):4d}" for v in pfm[b])
            lines.append(f"{b}  [ {row_str} ]")
        return "\n".join(lines) + "\n"
    if fmt_clean == "meme":
        calc = compute_ppm_and_pwm(pfm)
        ppm = calc["ppm"]
        n_sites = int(round(sum(calc["column_totals"]) / max(length, 1)))
        lines = [
            "MEME version 5",
            "",
            "ALPHABET= ACGT",
            "",
            "strands: + -",
            "",
            "Background letter frequencies",
            "A 0.25000 C 0.25000 G 0.25000 T 0.25000",
            "",
            f"MOTIF {matrix_id} {name}",
            f"letter-probability matrix: alength= 4 w= {length} nsites= {n_sites} E= 0",
        ]
        for i in range(length):
            lines.append(
                f"  {ppm['A'][i]:.6f}  {ppm['C'][i]:.6f}  {ppm['G'][i]:.6f}  {ppm['T'][i]:.6f}"
            )
        return "\n".join(lines) + "\n"
    if fmt_clean == "transfac":
        lines = [
            f"AC  {matrix_id}",
            "XX",
            f"ID  {name}",
            "XX",
            f"DE  {matrix_id} {name}; from JASPAR",
            "P0      A      C      G      T",
        ]
        for i in range(length):
            a_v = int(round(pfm["A"][i]))
            c_v = int(round(pfm["C"][i]))
            g_v = int(round(pfm["G"][i]))
            t_v = int(round(pfm["T"][i]))
            lines.append(f"{i + 1:02d}  {a_v:5d}  {c_v:5d}  {g_v:5d}  {t_v:5d}")
        lines.extend(["XX", "//"])
        return "\n".join(lines) + "\n"
    if fmt_clean == "yaml":
        lines = [
            f"matrix_id: {matrix_id}",
            f"name: {name}",
            f"collection: {data.get('collection', 'CORE')}",
            f"tax_group: {data.get('tax_group', 'vertebrates')}",
            "pfm:",
        ]
        for b in ("A", "C", "G", "T"):
            vals = ", ".join(str(int(round(v))) for v in pfm[b])
            lines.append(f"  {b}: [{vals}]")
        return "\n".join(lines) + "\n"

    raise ValueError(f"Unsupported format '{fmt}'. Choose from: json, jsonp, jaspar, meme, transfac, pfm, yaml.")


def emit_output(content: str, output_path: Optional[str]) -> None:
    """Writes content to file if output_path is specified; otherwise writes to stdout."""
    if output_path:
        target = Path(output_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")
    else:
        if len(content) > MAX_STDOUT_CHARS:
            content = content[:MAX_STDOUT_CHARS] + "\n... [Output truncated at 50,000 characters]"
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")


def cmd_resolve_tf_id(args: argparse.Namespace) -> Dict[str, Any]:
    """Resolves gene name and NCBI tax_id to JASPAR Matrix ID(s)."""
    verify_workspace_license(quiet=True)
    name = args.name.strip()
    tax_id = int(args.tax_id)

    if args.mock:
        fixture = load_mock_fixture(name)
        result = {
            "query": {"name": name, "tax_id": tax_id},
            "matrix_id": fixture["matrix_id"],
            "name": fixture["name"],
            "collection": fixture.get("collection", "CORE"),
            "count": 1,
            "results": [
                {
                    "matrix_id": fixture["matrix_id"],
                    "name": fixture["name"],
                    "collection": fixture.get("collection", "CORE"),
                    "base_id": fixture.get("base_id"),
                    "version": fixture.get("version"),
                    "species": fixture.get("species", [{"tax_id": tax_id, "name": "Homo sapiens"}]),
                }
            ],
            "citations": ACADEMIC_CITATION_NOTICE,
            "disclaimer": SCIENTIFIC_DISCLAIMER,
        }
    else:
        client = PoliteHttpClient(qps=10.0)
        query = urllib.parse.urlencode({"name": name, "tax_id": tax_id, "version": "latest"})
        api_resp = client.get_json(f"{JASPAR_API_BASE}/matrix/?{query}")
        results_list = api_resp.get("results", [])
        top_id = results_list[0]["matrix_id"] if results_list else None
        result = {
            "query": {"name": name, "tax_id": tax_id},
            "matrix_id": top_id,
            "name": name,
            "count": api_resp.get("count", len(results_list)),
            "results": results_list,
            "citations": ACADEMIC_CITATION_NOTICE,
            "disclaimer": SCIENTIFIC_DISCLAIMER,
        }

    emit_output(json.dumps(result, indent=2), args.output)
    return result


def cmd_get_tf_motif(args: argparse.Namespace) -> Dict[str, Any]:
    """Retrieves raw PFM counts and motif representation for a matrix ID."""
    verify_workspace_license(quiet=True)
    matrix_id = args.matrix_id.strip()

    if args.mock:
        data = load_mock_fixture(matrix_id)
    else:
        client = PoliteHttpClient(qps=10.0)
        data = client.get_json(f"{JASPAR_API_BASE}/matrix/{urllib.parse.quote(matrix_id)}/")

    data["citations"] = ACADEMIC_CITATION_NOTICE
    data["disclaimer"] = SCIENTIFIC_DISCLAIMER
    formatted = format_motif_output(data, args.format)
    emit_output(formatted, args.output)
    return data


def cmd_get_tf_metadata(args: argparse.Namespace) -> Dict[str, Any]:
    """Retrieves TF structural classification, species, and database cross-references."""
    verify_workspace_license(quiet=True)
    matrix_id = args.matrix_id.strip()

    if args.mock:
        raw = load_mock_fixture(matrix_id)
    else:
        client = PoliteHttpClient(qps=10.0)
        raw = client.get_json(f"{JASPAR_API_BASE}/matrix/{urllib.parse.quote(matrix_id)}/")

    metadata = {
        "matrix_id": raw.get("matrix_id"),
        "base_id": raw.get("base_id"),
        "version": raw.get("version"),
        "name": raw.get("name"),
        "collection": raw.get("collection", "CORE"),
        "tax_group": raw.get("tax_group", "vertebrates"),
        "species": raw.get("species", []),
        "superclass": raw.get("superclass", "Zinc-coordinating DNA-binding domains"),
        "class": raw.get("class", []),
        "family": raw.get("family", []),
        "subfamily": raw.get("subfamily", ""),
        "type": raw.get("type", "ChIP-seq"),
        "uniprot_ids": raw.get("uniprot_ids", []),
        "pdb_ids": raw.get("pdb_ids", []),
        "pubmed_ids": raw.get("pubmed_ids", []),
        "tffm": raw.get("tffm", {}),
        "citations": ACADEMIC_CITATION_NOTICE,
        "disclaimer": SCIENTIFIC_DISCLAIMER,
    }
    if args.format.lower() == "yaml":
        lines = [f"{k}: {json.dumps(v)}" for k, v in metadata.items()]
        emit_output("\n".join(lines) + "\n", args.output)
    else:
        emit_output(json.dumps(metadata, indent=2), args.output)
    return metadata


def cmd_get_tf_pwm(args: argparse.Namespace) -> Dict[str, Any]:
    """Fetches PFM and converts to PPM and log2-odds PWM."""
    verify_workspace_license(quiet=True)
    matrix_id = args.matrix_id.strip()

    if args.mock:
        raw = load_mock_fixture(matrix_id)
    else:
        client = PoliteHttpClient(qps=10.0)
        raw = client.get_json(f"{JASPAR_API_BASE}/matrix/{urllib.parse.quote(matrix_id)}/")

    calc = compute_ppm_and_pwm(
        raw["pfm"],
        pseudocount=args.pseudocount,
        gc_content=args.gc_content,
    )
    payload = {
        "matrix_id": raw.get("matrix_id"),
        "name": raw.get("name"),
        "length": calc["length"],
        "pseudocount": calc["pseudocount"],
        "gc_content": calc["gc_content"],
        "background": calc["background"],
        "pfm": raw["pfm"],
        "ppm": calc["ppm"],
        "pwm": calc["pwm"],
        "citations": ACADEMIC_CITATION_NOTICE,
        "disclaimer": SCIENTIFIC_DISCLAIMER,
    }
    emit_output(json.dumps(payload, indent=2), args.output)
    return payload


def cmd_infer_from_sequence(args: argparse.Namespace) -> Dict[str, Any]:
    """Matches protein amino acid sequence to candidate JASPAR matrices."""
    verify_workspace_license(quiet=True)
    sequence = args.sequence.strip().upper()
    if not sequence or len(sequence) < 5:
        raise ValueError("Protein sequence must be at least 5 amino acid residues.")

    if args.mock:
        payload = {
            "query_sequence_length": len(sequence),
            "count": 2,
            "results": [
                {
                    "matrix_id": "MA0139.1",
                    "name": "CTCF",
                    "e_value": 1.4e-42,
                    "identity_pct": 98.2,
                    "dbd_class": "Zinc finger, C2H2-related",
                },
                {
                    "matrix_id": "MA0488.2",
                    "name": "JUN",
                    "e_value": 3.8e-19,
                    "identity_pct": 84.5,
                    "dbd_class": "Basic leucine zipper factors (bZIP)",
                },
            ],
            "citations": ACADEMIC_CITATION_NOTICE,
            "disclaimer": SCIENTIFIC_DISCLAIMER,
        }
    else:
        client = PoliteHttpClient(qps=10.0)
        payload = client.get_json(f"{JASPAR_API_BASE}/infer/{urllib.parse.quote(sequence)}/")
        payload["citations"] = ACADEMIC_CITATION_NOTICE
        payload["disclaimer"] = SCIENTIFIC_DISCLAIMER

    emit_output(json.dumps(payload, indent=2), args.output)
    return payload


def cmd_get_tffm(args: argparse.Namespace) -> Dict[str, Any]:
    """Retrieves first-order and detailed HMM metadata for flexible models (TFFM)."""
    verify_workspace_license(quiet=True)
    tffm_id = args.tffm_id.strip()

    if args.mock:
        payload = {
            "tffm_id": tffm_id,
            "matrix_id": "MA0139.1",
            "name": "CTCF",
            "tax_group": "vertebrates",
            "log_p_1st_order": -198.4,
            "log_p_detailed": -212.7,
            "experiment_name": "ENCODE_ChIPseq_GM12878_CTCF",
            "citations": ACADEMIC_CITATION_NOTICE,
            "disclaimer": SCIENTIFIC_DISCLAIMER,
        }
    else:
        client = PoliteHttpClient(qps=10.0)
        payload = client.get_json(f"{JASPAR_API_BASE}/tffm/{urllib.parse.quote(tffm_id)}/")
        payload["citations"] = ACADEMIC_CITATION_NOTICE
        payload["disclaimer"] = SCIENTIFIC_DISCLAIMER

    emit_output(json.dumps(payload, indent=2), args.output)
    return payload


def build_parser() -> argparse.ArgumentParser:
    """Constructs the CLI argument parser supporting all six subcommands."""
    parser = argparse.ArgumentParser(
        prog="jaspar_api.py",
        description="JASPAR REST API v1 Client with Legal Guard, File Export, and Offline Mock Mode.",
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # 1. resolve_tf_id
    p_resolve = subparsers.add_parser("resolve_tf_id", help="Map gene name and tax_id to JASPAR Matrix ID.")
    p_resolve.add_argument("--name", required=True, help="Transcription factor gene symbol (e.g., CTCF, JUN).")
    p_resolve.add_argument("--tax-id", "--tax_id", dest="tax_id", type=int, default=9606, help="NCBI Taxonomy ID (default: 9606).")
    p_resolve.add_argument("--output", help="Optional file path to save JSON output.")
    p_resolve.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_resolve.set_defaults(func=cmd_resolve_tf_id)

    # 2. get_tf_motif
    p_motif = subparsers.add_parser("get_tf_motif", help="Retrieve PFM counts and motif representation.")
    p_motif.add_argument("--matrix-id", "--matrix_id", dest="matrix_id", required=True, help="JASPAR Matrix ID (e.g., MA0139.1).")
    p_motif.add_argument(
        "--format",
        default="json",
        choices=["json", "jsonp", "jaspar", "meme", "transfac", "pfm", "yaml"],
        help="Output format (default: json).",
    )
    p_motif.add_argument("--output", help="Optional file path to save output.")
    p_motif.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_motif.set_defaults(func=cmd_get_tf_motif)

    # 3. get_tf_metadata
    p_meta = subparsers.add_parser("get_tf_metadata", help="Retrieve TF class, family, species, and database links.")
    p_meta.add_argument("--matrix-id", "--matrix_id", dest="matrix_id", required=True, help="JASPAR Matrix ID (e.g., MA0139.1).")
    p_meta.add_argument("--format", default="json", choices=["json", "yaml"], help="Output format (default: json).")
    p_meta.add_argument("--output", help="Optional file path to save output.")
    p_meta.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_meta.set_defaults(func=cmd_get_tf_metadata)

    # 4. get_tf_pwm
    p_pwm = subparsers.add_parser("get_tf_pwm", aliases=["calculate_pwm"], help="Compute PPM and log2-odds PWM from PFM counts.")
    p_pwm.add_argument("--matrix-id", "--matrix_id", dest="matrix_id", required=True, help="JASPAR Matrix ID (e.g., MA0139.1).")
    p_pwm.add_argument("--pseudocount", type=float, default=0.8, help="Per-base pseudocount (default: 0.8).")
    p_pwm.add_argument("--gc-content", "--gc_content", dest="gc_content", type=float, default=0.50, help="Background GC content (default: 0.50).")
    p_pwm.add_argument("--output", help="Optional file path to save output.")
    p_pwm.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_pwm.set_defaults(func=cmd_get_tf_pwm)

    # 5. infer_from_sequence
    p_infer = subparsers.add_parser("infer_from_sequence", aliases=["infer_motifs"], help="Match protein sequence to candidate matrices.")
    p_infer.add_argument("--sequence", required=True, help="Amino acid sequence string.")
    p_infer.add_argument("--tax-id", "--tax_id", dest="tax_id", type=int, default=9606, help="NCBI Taxonomy ID (default: 9606).")
    p_infer.add_argument("--output", help="Optional file path to save output.")
    p_infer.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_infer.set_defaults(func=cmd_infer_from_sequence)

    # 6. get_tffm
    p_tffm = subparsers.add_parser("get_tffm", aliases=["get_tffm_profile"], help="Retrieve first-order HMM TFFM metadata.")
    p_tffm.add_argument("--tffm-id", "--tffm_id", dest="tffm_id", required=True, help="JASPAR TFFM ID (e.g., TFFM0089.3).")
    p_tffm.add_argument("--output", help="Optional file path to save output.")
    p_tffm.add_argument("--mock", action="store_true", help="Use deterministic offline mock fixtures.")
    p_tffm.set_defaults(func=cmd_get_tffm)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entrypoint for jaspar_api.py."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        args.func(args)
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        sys.stderr.write(f"[VALIDATION ERROR] {exc}\n")
        return 2
    except Exception as exc:
        sys.stderr.write(f"[EXECUTION ERROR] {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
