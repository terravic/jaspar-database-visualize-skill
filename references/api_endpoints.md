# JASPAR REST API v1 Endpoint Reference

## Base URL and Rate Limiting
- **Base Endpoint**: `https://jaspar.elixir.no/api/v1/`
- **Rate Limit Policy**: Maximum 10 queries per second (10 QPS) with exponential backoff on HTTP 429 and 5xx status codes.
- **License Terms**: Creative Commons Attribution 4.0 International (CC BY 4.0).

---

## 1. Matrix Listing and TF Resolution (`GET /api/v1/matrix/`)
Resolves transcription factor names, species taxonomy identifiers, and structural families to canonical JASPAR Matrix IDs.

### Query Parameters
| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `name` | string | Case-insensitive transcription factor symbol | `CTCF`, `JUN`, `TP53` |
| `tax_id` | integer | NCBI Taxonomy Identifier | `9606` (Homo sapiens), `10090` (Mus musculus) |
| `collection` | string | JASPAR collection filter | `CORE`, `UNVALIDATED` |
| `tax_group` | string | Taxonomic supergroup filter | `vertebrates`, `plants`, `insects`, `fungi` |
| `version` | string | Return only the latest profile version | `latest` |
| `page_size` | integer | Number of records per page | `25` |

### Sample Request
```text
GET https://jaspar.elixir.no/api/v1/matrix/?name=CTCF&tax_id=9606&version=latest
```

---

## 2. Matrix Profile and PFM Retrieval (`GET /api/v1/matrix/{matrix_id}/`)
Retrieves the complete Position Frequency Matrix (PFM), TFClass structural annotations, UniProt links, and experimental evidence for a specific matrix ID.

### Supported Content Types and Formats
- `json`: Structured JSON object containing `pfm` (`A`, `C`, `G`, `T` arrays), `class`, `family`, `species`, `uniprot_ids`, and `pubmed_ids`.
- `jaspar`: Header line (`>MA0139.1 CTCF`) followed by bracketed A, C, G, T count rows.
- `meme`: MEME suite minimal motif format (`MEME version 5`, background frequencies, and `letter-probability matrix`).
- `transfac`: TRANSFAC table format (`AC`, `ID`, `DE`, `P0`, and positional base counts).
- `yaml`: YAML serialized representation of the matrix and metadata.
- `pfm`: Four lines of whitespace-separated nucleotide counts ordered A, C, G, T.

---

## 3. Amino Acid Homology Inference (`GET /api/v1/infer/{sequence}/`)
Infers candidate transcription factor binding profiles by aligning a query protein amino acid sequence against DNA-binding domain sequences of known JASPAR TFs.

### Path Parameter
- `sequence`: Raw single-letter amino acid sequence string (minimum 20 residues recommended).

---

## 4. Transcription Factor Flexible Models (`GET /api/v1/tffm/{tffm_id}/`)
Retrieves first-order and detailed Hidden Markov Model (HMM) metadata capturing dinucleotide positional dependencies (`TFFM`).

### Path Parameter
- `tffm_id`: JASPAR TFFM identifier (e.g., `TFFM0001.1`).
