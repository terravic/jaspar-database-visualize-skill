# TFClass Structural and Taxonomic Classification Reference

JASPAR categorizes transcription factors according to the hierarchical **TFClass** structural framework (`Superclass -> Class -> Family -> Subfamily`), based on the three-dimensional fold and sequence homology of their DNA-binding domains (DBDs).

---

## 1. Major TFClass Superclasses Represented in JASPAR

| Superclass ID | Superclass Name | Representative Classes | Canonical Examples in Skill Fixtures |
| :--- | :--- | :--- | :--- |
| `1` | Basic domains | Basic leucine zipper factors (bZIP), Basic helix-loop-helix factors (bHLH) | `JUN` (`MA0488.2`, bZIP AP-1 complex) |
| `2` | Zinc-coordinating DNA-binding domains | C2H2 zinc finger factors, Nuclear receptors with C4 zinc fingers | `CTCF` (`MA0139.1`, C2H2 11-zinc finger insulator) |
| `3` | Helix-turn-helix domains | Homeo domain factors, Fork head / winged helix factors, Tryptophan cluster factors | `FOXA1`, `HOXA1`, `SPI1` |
| `4` | All-alpha-helical DNA-binding domains | High-mobility group (HMG) domain factors, MADS box factors | `SOX2` (`MA0143.4`, SRY-related HMG box) |
| `5` | Alpha-helices exposed by beta-structures | SAND domain factors, E2F/DP family | `E2F1` |
| `6` | Immunoglobulin fold | Rel homology region (RHR) factors, STAT domain factors, p53 domain factors | `TP53` (`MA0106.3`, p53 tetrameric DNA-binding domain) |
| `7` | Beta-Beta-Alpha-Me finger | GCM domain factors | `GCM1` |
| `8` | Beta-sheet binding to DNA | TATA-binding proteins (TBP) | `TBP` |

---

## 2. Reference Profiles Included in Deterministic Mock Mode

1. **`MA0139.1` (`CTCF` - CCCTC-Binding Factor)**:
   - **Superclass**: Zinc-coordinating DNA-binding domains
   - **Class**: Zinc finger, C2H2-related
   - **Family**: More than 3 adjacent zinc finger factors
   - **UniProt ID**: `P49711`
   - **Species**: *Homo sapiens* (`9606`)
   - **Biological Role**: Chromatin architecture protein mediating topologically associating domain (TAD) boundaries and enhancer-promoter insulation.

2. **`MA0488.2` (`JUN` - Transcription Factor Jun / AP-1)**:
   - **Superclass**: Basic domains
   - **Class**: Basic leucine zipper factors (bZIP)
   - **Family**: Jun-related factors
   - **UniProt ID**: `P05412`
   - **Species**: *Homo sapiens* (`9606`)
   - **Biological Role**: Binds palindromic 12-O-tetradecanoylphorbol-13-acetate (TPA) response elements (`TGAGTCA`) as homodimers or JUN/FOS heterodimers.

3. **`MA0106.3` (`TP53` - Cellular Tumor Antigen p53)**:
   - **Superclass**: Immunoglobulin fold
   - **Class**: p53 domain factors
   - **Family**: p53-related factors
   - **UniProt ID**: `P04637`
   - **Species**: *Homo sapiens* (`9606`)
   - **Biological Role**: Tumor suppressor binding symmetrical decameric half-sites (`RRRCWWGYYY`) to regulate cell cycle arrest and apoptosis.

4. **`MA0143.4` (`SOX2` - Transcription Factor SOX-2)**:
   - **Superclass**: All-alpha-helical DNA-binding domains
   - **Class**: High-mobility group (HMG) domain factors
   - **Family**: SOX-related factors
   - **UniProt ID**: `P48431`
   - **Species**: *Homo sapiens* (`9606`)
   - **Biological Role**: Core pluripotency transcription factor bending the DNA minor groove at `CTTTGTT` consensus elements.
