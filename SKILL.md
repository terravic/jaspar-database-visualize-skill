---
name: jaspar-database-visualize
description: Query the JASPAR REST API v1 for transcription factor binding profiles, compute PFM/PPM/PWM matrices, render vector SVG and HTML5 sequence logos, scan DNA sequences across forward and reverse strands, and compile standalone interactive HTML5 dashboards.
---

# JASPAR Database Visualize (`jaspar-database-visualize`)

## 1. Purpose & Scope

`jaspar-database-visualize` is a computational biology skill for querying, analyzing, and visualizing transcription factor (TF) DNA-binding motifs from the JASPAR open-access database (`https://jaspar.elixir.no/api/v1/`).

### Core Analytical Capabilities
1. **JASPAR REST API v1 Retrieval & Offline Fallback (`scripts/jaspar_api.py`)**:
   - Resolves human and model-organism TF gene symbols (`CTCF`, `TP53`, `MYC`, `FOXA1`, `SOX2`, `JUN`, `GATA1`) to versioned JASPAR Matrix IDs (`MA0139.1`, `MA0106.3`, `MA0147.3`, `MA0148.4`, `MA0143.4`, `MA0488.2`, `MA0035.4`).
   - Retrieves raw Position Frequency Matrices (PFMs), TFClass structural domain hierarchies, UniProt accessions, PDB structural identifiers, and species taxonomy.
   - Computes Laplace/Dirichlet pseudocount-normalized Position Probability Matrices (PPMs) and background-adjusted $\log_2$-odds Position Weight Matrices (PWMs).
   - Exports motifs in `JSON`, `JASPAR`, `MEME`, `TRANSFAC`, `PFM`, and `YAML` formats.
2. **Information-Theoretic Vector Sequence Logos (`scripts/visualize_logo.py`)**:
   - Computes per-position Shannon entropy ($H_i$), Schneider & Stephens small-sample error correction ($e_n$), and information content ($R_i = 2 - (H_i + e_n)$ bits).
   - Renders standalone vector `.svg` files and interactive HTML5 widgets across Forward ($5' \to 3'$) and Reverse Complement ($3' \leftarrow 5'$) orientations, with IUPAC degenerate consensus string derivation (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`).
3. **PFM / PPM / PWM Matrix Heatmap Explorer (`scripts/visualize_matrix.py`)**:
   - Generates color-scaled heatmaps and interactive HTML5 widgets with live client-side sliders for pseudocount ($p \in [0.1, 5.0]$) and background GC fraction ($GC \in [0.20, 0.80]$).
4. **Dual-Strand DNA Sequence Scanner (`scripts/visualize_scan.py`)**:
   - Scans user-provided FASTA files or raw DNA strings across both Forward (`+`) and Reverse (`-`) strands using sliding-window PWM scoring, normalizing raw log-odds scores to $[0.00, 1.00]$ relative affinity (`Optimal`, `High`, `Moderate`) and rendering a genomic coordinate locus map.
5. **Standalone Interactive HTML5 Dashboard (`scripts/generate_dashboard.py` & `scripts/run_pipeline.py`)**:
   - Compiles processed motif metadata, vector sequence logos, PFM/PPM/PWM heatmaps, structural TFClass hierarchy, familial cluster comparisons, dual-strand sequence scan hits, and 3D interactive views into a single self-contained `dashboard.html` file.

---

## 2. Operational Execution Flow

### Step 1: Process Input Data via `scripts/run_pipeline.py`
Process the target TF, custom motif JSON, or DNA sequence (`--sequence` or `--fasta`) first so that all downstream outputs and `dashboard.html` reflect the input dataset:

```bash
# Run full pipeline on a specified TF and inline DNA sequence
python3 scripts/run_pipeline.py \
  --name "MYC" \
  --tax-id 9606 \
  --sequence "GCTAGCTAGCGACCACGTGGTCGATCGATCGGACCACGTGGTCGATCGATCG" \
  --threshold 0.80 \
  --output-dir outputs/myc_pipeline_output

# Run with --mock for deterministic offline execution
python3 scripts/run_pipeline.py \
  --name "CTCF" \
  --tax-id 9606 \
  --fasta sample_data/sample_enhancer.fasta \
  --threshold 0.80 \
  --output-dir outputs/ctcf_pipeline_output \
  --mock
```

### Step 2: Module CLI Invocations

1. **`scripts/jaspar_api.py`** -- Resolve TF IDs, fetch PFM/metadata, compute PWM, or infer motifs:
   ```bash
   # Resolve TF symbol to latest JASPAR Matrix ID
   python3 scripts/jaspar_api.py resolve_tf_id \
     --name "CTCF" --tax-id 9606 --output outputs/ctcf_resolve.json --mock

   # Fetch PFM matrix in JSON, JASPAR, MEME, TRANSFAC, PFM, or YAML format
   python3 scripts/jaspar_api.py get_tf_motif \
     --matrix-id "MA0139.1" --format json --output outputs/ctcf_motif.json --mock

   # Fetch TFClass structural and taxonomic metadata
   python3 scripts/jaspar_api.py get_tf_metadata \
     --matrix-id "MA0139.1" --output outputs/ctcf_metadata.json --mock

   # Compute PPM and log2-odds PWM with custom pseudocount and GC fraction
   python3 scripts/jaspar_api.py get_tf_pwm \
     --matrix-id "MA0139.1" --pseudocount 0.8 --gc-content 0.50 --output outputs/ctcf_pwm.json --mock

   # Infer candidate JASPAR motifs from a protein sequence
   python3 scripts/jaspar_api.py infer_from_sequence \
     --sequence "YACPVESCDRRFSDSSNLTRHIRIHTGQKP" --tax-id 9606 \
     --output outputs/infer_hits.json --mock

   # Query TFFM first-order Markov dinucleotide models
   python3 scripts/jaspar_api.py get_tffm \
     --tffm-id "TFFM0089.3" --output outputs/tffm_meta.json --mock
   ```

2. **`scripts/visualize_logo.py`** -- Render vector SVG and HTML5 sequence logos:
   ```bash
   python3 scripts/visualize_logo.py \
     --matrix-json outputs/ctcf_motif.json \
     --output-svg outputs/ctcf_logo.svg \
     --output-html outputs/ctcf_logo.html \
     --mock
   ```

3. **`scripts/visualize_matrix.py`** -- Render PFM/PPM/PWM matrix heatmaps and exports:
   ```bash
   python3 scripts/visualize_matrix.py \
     --matrix-json outputs/ctcf_motif.json \
     --pseudocount 0.8 \
     --gc-content 0.50 \
     --output-json outputs/ctcf_matrix.json \
     --output-html outputs/ctcf_matrix.html \
     --mock
   ```

4. **`scripts/visualize_scan.py`** -- Scan FASTA or raw DNA sequence across `+` and `-` strands:
   ```bash
   python3 scripts/visualize_scan.py \
     --matrix-json outputs/ctcf_motif.json \
     --fasta sample_data/sample_enhancer.fasta \
     --threshold 0.80 \
     --output-json outputs/ctcf_hits.json \
     --output-html outputs/ctcf_hits.html \
     --mock
   ```

5. **`scripts/generate_dashboard.py`** -- Compile standalone interactive `dashboard.html`:
   ```bash
   python3 scripts/generate_dashboard.py \
     --matrix-id "MA0139.1" \
     --motif-json outputs/ctcf_motif.json \
     --scan-json outputs/ctcf_hits.json \
     --output outputs/ctcf_dashboard.html \
     --mock
   ```

---

## 3. Interactive HTML5 Dashboard Protocol

When generating outputs for a query, compile and save `dashboard.html` reflecting the processed input data:

1. **Generate the Dashboard with Processed Data**: Run `scripts/run_pipeline.py` or `scripts/generate_dashboard.py` with the `--motif-json`, `--scan-json`, `--fasta`, or `--sequence` parameters.
2. **Save Local Workspace Copy**: Write the compiled dashboard to `./dashboard.html` in the project root directory.
3. **Interactive Components Included in `dashboard.html`**:
   - **Light / Dark / System Theme Selector**: Top-right segmented control (`#theme-btn-light`, `#theme-btn-dark`, `#theme-btn-system`) switching between light, dark, and automatic OS appearance modes.
   - **Input Data & Live Sequence Scanner**: Displays processed sequence metrics (length, GC%, base composition, total hits) and provides an editable sequence input with a **Run Sequence Scan** button that recalculates forward/reverse strand hits in real time.
   - **Interactive Relationship Node Map & Structural Data Diagram**: Interactive 3D SVG viewer (`#viewport3DSurface`) supporting mouse rotation, Shift+drag pan, scroll zoom, and node inspection for both TF cluster topology and 3D B-DNA double-helix groove geometry.
   - **Sequence Logo & PFM/PPM/PWM Matrix Explorer**: Live toggles for forward/reverse-complement orientation, bits/probability scaling, pseudocount/GC sliders, and `JASPAR`/`MEME`/`TRANSFAC`/`YAML` exports.

---

## 4. Attribution & Scientific Disclaimer

Every output and dashboard generated by this skill includes:
1. **JASPAR 2026 & 2024 Citations**:
   - Ovek Baydar, D. et al. *JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models.* Nucleic Acids Research 54, D184-D193 (2026), `doi:10.1093/nar/gkaf1209`.
   - Rauluseviciute, I. et al. *JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles.* Nucleic Acids Research 52, D174-D182 (2024), `doi:10.1093/nar/gkad1059`.
2. **License Attribution**: Creative Commons Attribution 4.0 International (`CC BY 4.0`).
3. **Scientific Disclaimer**: *"Predicted binding sites must be experimentally validated via ChIP-seq, CUT&RUN, or reporter assays."*
