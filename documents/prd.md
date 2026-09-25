# Product Requirement Document (PRD): JASPAR Database Visualize (`jaspar-database-visualize`)

## 1. Document Metadata & Overview
- **Skill Identifier**: `jaspar-database-visualize`
- **Version**: `2.2.0`
- **Domain**: Computational Biology / Gene Regulation / Transcription Factor Genomics
- **Primary Upstream Data Source**: JASPAR REST API v1 (`https://jaspar.elixir.no/api/v1/`)
- **Data License**: Creative Commons Attribution 4.0 International (`CC BY 4.0`)
- **Software License**: Apache License, Version 2.0 (`LICENSE`)

---

## 2. Target Audience & User Groups

The `jaspar-database-visualize` toolkit supports two primary user groups:

1. **Biomedical Researchers, Clinicians, and Students**:
   - **Workflow**: Supply a transcription factor gene symbol (such as `MYC`, `TP53`, `CTCF`, `FOXA1`, or `SOX2`) and optionally a raw DNA sequence or `.fasta` file from a laboratory experiment or variant report.
   - **Core Need**: Obtain a self-contained interactive visual report (`dashboard.html`) that highlights where a transcription factor binds along their DNA sequence, displays the structural family of the protein, and enables interactive testing of sequence variants in any web browser.
2. **Computational Biologists and Bioinformatics Engineers**:
   - **Workflow**: Execute the modular CLI scripts (`scripts/run_pipeline.py`, `scripts/jaspar_api.py`, `scripts/visualize_logo.py`, `scripts/visualize_matrix.py`, `scripts/visualize_scan.py`, `scripts/generate_dashboard.py`) in genomics pipelines or batch analyses.
   - **Core Need**: Verified Shannon information content ($R_i$) calculations, configurable pseudocount/GC-background Position Weight Matrices (PWMs), dual-strand sliding-window sequence scanning, and standardized motif format exports (`JASPAR`, `MEME`, `TRANSFAC`, `PFM`, `YAML`, `JSON`).

---

## 3. Repository Structure & File Organization

All file paths in the repository are organized relative to the project root directory:

- `LICENSE`: Apache License, Version 2.0 text.
- `SKILL.md`: Operational specification and CLI reference for the skill.
- `SKILL_LICENSES.md` & `.licenses/jaspar_database_visualize_LICENSE.txt`: Upstream JASPAR `CC BY 4.0` data license attribution records.
- `README.md`: Documentation, non-technical user guide with real-world examples, relative directory tree, and CLI usage guide.
- `dashboard.html`: Pre-compiled standalone interactive HTML5 dashboard at the project root.
- `pytest.ini` & `requirements.txt`: Test runner configuration (`--import-mode=importlib -p no:cacheprovider`) and Python package dependencies.
- `assets/`:
  - `dashboard_template.html`: Baseline interactive HTML5 dashboard template.
  - `jaspar_skill_workflow.png`: Technical architecture and data-flow diagram.
  - `motif_matrix_schema.json` & `scan_results_schema.json`: JSON Schema (Draft 2020-12) validation definitions.
- `documents/prd.md`: Product Requirement Document (this specification).
- `references/`:
  - `api_endpoints.md`: JASPAR REST API v1 endpoint and rate-limit documentation.
  - `citation.bib`: BibTeX academic references for JASPAR 2026 and JASPAR 2024.
  - `dashboard_schema.md`: Schema and UI specification for `window.JASPAR_DASHBOARD_STATE` in `dashboard.html`.
  - `tf_classification.md`: TFClass structural superclass, class, and family reference.
- `sample_data/`:
  - `ctcf_ma0139_1.json`, `jun_ma0488_2.json`, `tp53_ma0106_3.json`, `sox2_ma0143_4.json`: Reference PFM and metadata fixtures.
  - `sample_enhancer.fasta`: 1,000 bp synthetic human enhancer sequence containing known forward and reverse binding sites.
- `scripts/`:
  - `jaspar_api.py`, `visualize_logo.py`, `visualize_matrix.py`, `visualize_scan.py`, `generate_dashboard.py`, `run_pipeline.py`.
- `tests/`:
  - `run_tests.py` & `test_jaspar_visualize.py`: Automated verification suite.

---

## 4. Scope & Functional Modules

### 4.1 Module 1: JASPAR REST API v1 Client & Offline Fallback (`scripts/jaspar_api.py`)
- **TF Symbol Resolution (`resolve_tf_id`)**: Resolves gene symbols (`CTCF`, `TP53`, `MYC`, `FOXA1`, `SOX2`, `JUN`, `GATA1`, etc.) and NCBI taxonomy filters (`--tax-id 9606` for *Homo sapiens*) to the latest versioned JASPAR Matrix ID (`MA0139.1`, `MA0106.3`, `MA0147.3`, `MA0148.4`, `MA0143.4`, `MA0488.2`, `MA0035.4`).
- **Motif & Metadata Retrieval (`get_tf_motif`, `get_tf_metadata`)**: Fetches raw Position Frequency Matrices (PFMs), TFClass structural hierarchy (`Superclass -> Class -> Family -> Subfamily`), UniProt accessions, and PDB structural identifiers.
- **PWM Calculation (`get_tf_pwm` / `calculate_pwm`)**: Converts raw PFM counts into Laplace/Dirichlet pseudocount-adjusted Position Probability Matrices (PPMs) and $\log_2$-odds Position Weight Matrices (PWMs) against user-specified background GC fractions:
  $$P_{b,i} = \frac{C_{b,i} + p}{N_i + 4p}, \quad W_{b,i} = \log_2\left(\frac{P_{b,i}}{B_b}\right)$$
- **Sequence-to-Motif Inference & TFFM Support (`infer_from_sequence` / `infer_motifs`, `get_tffm` / `get_tffm_profile`)**: Queries candidate motifs matching protein sequences and retrieves Transcription Factor Flexible Model (TFFM) first-order Markov dinucleotide profiles.
- **Deterministic Offline Mode (`--mock`)**: Provides local fallback fixtures in `sample_data/` (`CTCF`, `JUN`, `TP53`, `SOX2`, `sample_enhancer.fasta`) plus built-in canonical profiles (`MYC`, `FOXA1`, `GATA1`) so all pipelines and tests execute deterministically without network access.

### 4.2 Module 2: Vector Sequence Logo Engine (`scripts/visualize_logo.py`)
- Computes per-column total counts $N_i = \sum_{b \in \{A,C,G,T\}} C_{b,i}$, Shannon entropy $H_i = -\sum_{b} P_{b,i} \log_2(P_{b,i})$, small-sample correction $e_n = \frac{3}{2 \ln(2) N_i}$, and information content $R_i = \max(0, 2 - (H_i + e_n))$ in bits ($[0, 2]$).
- Stacks vector `<path>` nucleotide glyphs (`A`: `#10B981`, `C`: `#3B82F6`, `G`: `#F59E0B`, `T`: `#EF4444`) sorted ascending by height $h_{b,i} = P_{b,i} \cdot R_i$.
- Supports both Forward ($5' \to 3'$) and Reverse Complement ($3' \leftarrow 5'$) orientations and derives the 15-letter IUPAC degenerate consensus string (`A`, `C`, `G`, `T`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`).

### 4.3 Module 3: Matrix Heatmap & Format Exporter (`scripts/visualize_matrix.py`)
- Computes PFM, PPM, and PWM matrices along with per-column entropy and information content vectors.
- Exports formatted motif representations in `JASPAR`, `MEME` (v5), `TRANSFAC`, `PFM`, `YAML`, and `JSON`.
- Renders an interactive HTML5 heatmap widget with client-side sliders for pseudocount ($0.1 \text{ to } 5.0$) and background GC content ($0.20 \text{ to } 0.80$).

### 4.4 Module 4: Dual-Strand DNA Sequence Scanner (`scripts/visualize_scan.py`)
- Accepts user FASTA files (`--fasta`), raw inline DNA sequence strings (`--sequence`), or custom motif JSON files (`--motif-json`) outside `sample_data/`.
- Evaluates every sliding window of length $L$ across both Forward (`+`) and Reverse (`-`) strands using the $\log_2$-odds PWM.
- Normalizes raw scores $S$ to relative affinity scores:
  $$\text{Relative Score} = \frac{S - S_{\min}}{S_{\max} - S_{\min}} \in [0.0, 1.0]$$
- Categorizes hits above the user threshold (`--threshold`, default `0.80`) into `Optimal` ($\ge 0.95$), `High` ($0.85 - 0.95$), and `Moderate` ($0.80 - 0.85$) tiers, and preserves the cleaned input sequence (`scan_data["sequence"]`) for downstream dashboard compilation.

### 4.5 Module 5: Unified Pipeline & Interactive HTML5 Dashboard (`scripts/run_pipeline.py` & `scripts/generate_dashboard.py`)
- Compiles processed motif metadata, vector SVG sequence logos, PFM/PPM/PWM matrices, TFClass classifications, familial cluster comparisons, and dual-strand sequence scan results into a standalone `dashboard.html` file.
- Embeds the processed state into `window.JASPAR_DASHBOARD_STATE`, ensuring that `dashboard.html` reflects the user-supplied TF profile and DNA sequence.
- Interactive Dashboard Components:
  - **Light / Dark Theme Switcher (`#themeToggleBtn`)**: Text-based header control (`Theme: Light` / `Theme: Dark`) that switches the dashboard palette between Light and Dark modes.
  - **Input Data & Live Sequence Scanner**: Displays sequence statistics (length in bp, observed GC%, A/C/G/T base counts, total detected hits), an editable DNA sequence box, preset sequence loaders, a `+ Inject Motif Site` helper, a relative score threshold slider, and a **Run Sequence Scan** button with status feedback.
  - **3D Interactive SVG Viewport (`#viewport3DSurface`)**: Renders an interactive 3D SVG **Motif Relationship Network** (linking the active TF hub to its TFClass domain, UniProt/PDB records, familial cluster paralogs, and top genomic binding hits) and a **3D DNA Binding Groove** (rendering the 3D B-DNA double-helix base pairs scaled by per-position bit information content $R_i$) with mouse drag rotation, scroll zoom, and node click inspection.

---

## 5. End-to-End Architecture & Data Flow

```text
+-----------------------------------------------------------------------------+
|                                USER INPUTS                                  |
|  - Target TF Symbol (e.g., MYC, FOXA1, CTCF, TP53) or Matrix ID             |
|  - Custom User DNA Sequence (--sequence) or External FASTA File (--fasta)   |
|  - Optional Custom PFM Matrix JSON (--motif-json)                           |
+--------------------------------------+--------------------------------------+
                                       |
                                       v
+-----------------------------------------------------------------------------+
|             STAGE 1: DATA INGESTION & MATHEMATICAL PROCESSING               |
|  1. scripts/jaspar_api.py      -> Resolves Matrix ID, fetches PFM/TFClass   |
|  2. scripts/visualize_logo.py  -> Computes H_i, e_n, R_i, IUPAC, SVG logos  |
|  3. scripts/visualize_matrix.py-> Computes PPM, PWM, MEME/TRANSFAC exports  |
|  4. scripts/visualize_scan.py  -> Scans (+) and (-) strands, records hits   |
+--------------------------------------+--------------------------------------+
                                       |
                                       v
+-----------------------------------------------------------------------------+
|             STAGE 2: DYNAMIC HTML5 DASHBOARD COMPILATION                    |
|  - scripts/generate_dashboard.py embeds window.JASPAR_DASHBOARD_STATE       |
|  - Compiles standalone interactive HTML5 dashboard (dashboard.html)         |
|  - Saves output to outputs/<run>/dashboard.html and ./dashboard.html        |
+-----------------------------------------------------------------------------+
```

---

## 6. Edge Cases & Input Validation Rules

1. **Sequence Length Boundaries**:
   - Sequences shorter than the motif length ($L_{\text{seq}} < L_{\text{motif}}$) return `0` hits with an explicit diagnostic warning in `warnings` rather than raising an unhandled exception.
   - Sequences exceeding `100,000 bp` (`100 kb`) raise a descriptive `ValueError` instructing the user to segment the locus into $\le 100\text{ kb}$ windows.
2. **Ambiguous Base Handling (`N`) in Target DNA**:
   - Ambiguous `N` bases in scanned DNA sequences are assigned a neutral log-odds weight ($0.0\text{ bits}$) at that column, while unsupported non-DNA characters raise a `ValueError`.
3. **Background GC Fraction Bounds**:
   - Background GC content is validated within $[0.01, 0.99]$ (`--gc-content`) and pseudocounts within $(0.0, 100.0]$ (`--pseudocount`) to prevent $\log_2(0)$ division errors during PWM derivation.
4. **License Verification**:
   - Every script verifies or initializes `.licenses/jaspar_database_visualize_LICENSE.txt` (`CC BY 4.0` JASPAR data attribution) and `SKILL_LICENSES.md` prior to execution.

---

## 7. Verification & Acceptance Criteria

All components are verified by the automated `pytest` suite (`tests/test_jaspar_visualize.py`):
1. `test_legal_compliance_and_apache_headers`: Verifies `LICENSE`, `SKILL_LICENSES.md`, `.licenses/jaspar_database_visualize_LICENSE.txt`, and Apache 2.0 source headers across all scripts.
2. `test_ppm_and_pwm_mathematics`: Verifies Laplace/Dirichlet pseudocount normalization and $\log_2$-odds PWM accuracy.
3. `test_shannon_entropy_and_iupac_consensus`: Verifies Shannon entropy $H_i$, small-sample correction $e_n$, bit heights $R_i$, reverse-complement symmetry, and IUPAC degenerate mapping.
4. `test_dual_strand_sequence_scanner_and_edge_cases`: Verifies forward (`+`) and reverse (`-`) strand hit detection, short-sequence warning handling, and the `>100 kb` sequence length check.
5. `test_dashboard_generation_and_mandatory_citations`: Verifies `dashboard.html` state injection, JASPAR 2026/2024 citations, CC BY 4.0 attribution, and the experimental validation disclaimer.
6. `test_cli_acceptance_suite`: Executes all CLI subcommands end-to-end.
7. `test_custom_user_data_pipeline_and_dashboard_reflection`: Verifies that custom user DNA sequences and TFs outside `sample_data/` (`MYC` / `MA0147.3`) are processed and embedded directly into `dashboard.html`.
