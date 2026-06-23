# c-Data_generated — Generated Data

This folder contains all data files produced by the project's analysis pipeline. Files are grouped below by processing stage, following the order in which they are generated. Each entry describes the producing script, inputs, and format.

---

## Stage 1 — PolyX detection

### `polyx.txt`

**Description:** Full catalogue of polyX regions detected across all filtered protein sequences. This is the central reference file for all polyX coordinates and amino acid types used in subsequent steps.

**Produced by:** `1-polyx2_standalone.pl` (see `d-Software/README.md`)

**Input:** `filtered_seqs.fasta` (`b-Data_processed`); parameters: minimum = 8, window = 10.

**Format:** Tab-separated, with header.

| Column | Description |
|---|---|
| `ID` | Protein identifier |
| `Start` | 1-based start coordinate of the polyX region |
| `End` | 1-based end coordinate of the polyX region |
| `Aa` | Amino acid forming the polyX run |
| `polyX` | Sequence of the polyX region |

**Number of entries:** 80,446 polyX regions.

---

## Stage 2 — Ortholog pair classification

### `11polyx.txt`

**Description:** Subset of filtered ortholog pairs where both proteins carry exactly one polyX region. Serves as the input for pairwise alignment and one-to-one classification.

**Produced by:** `4-split_pairs.pl`

**Input:** `filtered_pairs.txt` (`b-Data_processed`), `polyx.txt`.

**Format:** Tab-separated, no header.

| Column | Description |
|---|---|
| 1 | Protein ID 1 |
| 2 | Protein ID 2 |
| 3 | Normalized similarity score |
| 4 | PolyX count per protein, format `[count1\|count2]` |

**Number of pairs:** 89,845.

---

### `m1polyx.txt`

**Description:** Subset of filtered ortholog pairs where one protein has exactly one polyX and the other has more than one.

**Produced by:** `4-split_pairs.pl` (same run as `11polyx.txt`)

**Format:** Same four-column tab-separated format as `11polyx.txt`.

**Number of pairs:** 60,060.

---

### `mmpolyx.txt`

**Description:** Subset of filtered ortholog pairs where both proteins carry more than one polyX region.

**Produced by:** `4-split_pairs.pl` (same run as `11polyx.txt`)

**Format:** Same four-column tab-separated format as `11polyx.txt`.

**Number of pairs:** 25,694.

---

### `01polyx.txt`

**Description:** Subset of filtered ortholog pairs where one protein has exactly one polyX and the other has none.

**Produced by:** `4-split_pairs.pl` (same run as `11polyx.txt`)

**Format:** Same four-column tab-separated format as `11polyx.txt`.

**Number of pairs:** 811,796.

---

### `00polyx.txt`

**Description:** Subset of filtered ortholog pairs where neither protein carries any polyX region.

**Produced by:** `4-split_pairs.pl` (same run as `11polyx.txt`)

**Format:** Same four-column tab-separated format as `11polyx.txt`.

**Number of pairs:** 7,479,008.

---

### `m0polyx.txt`

**Description:** Subset of filtered ortholog pairs where one protein has more than one polyX and the other has none.

**Produced by:** `4-split_pairs.pl` (same run as `11polyx.txt`)

**Format:** Same four-column tab-separated format as `11polyx.txt`.

**Number of pairs:** 237,860.

---

### `polyx_pair.txt`

**Description:** Annotated version of `11polyx.txt` in which each pair is enriched with the polyX coordinates and amino acid type for both proteins, retrieved from `polyx.txt`.

**Produced by:** `5-get_info.pl`

**Input:** `11polyx.txt` (this folder), `polyx.txt` (this folder).

**Format:** Tab-separated, with header.

| Column | Description |
|---|---|
| `ID1` | Protein ID 1 |
| `Aa1` | Amino acid forming the polyX of protein 1 |
| `Start1` | PolyX start coordinate in protein 1 |
| `End1` | PolyX end coordinate in protein 1 |
| `ID2` | Protein ID 2 |
| `Aa2` | Amino acid forming the polyX of protein 2 |
| `Start2` | PolyX start coordinate in protein 2 |
| `End2` | PolyX end coordinate in protein 2 |

**Number of pairs:** 89,845.

---

## Stage 3 — Alignment-based overlap classification

### `desglose_onetoone.txt`

**Description:** Classification of all one-to-one ortholog pairs (from `polyx_pair.txt`) into four categories based on pairwise sequence alignment and polyX overlap. Each pair is labelled A (same amino acid, overlapping), B (same amino acid, not overlapping), C (different amino acids, overlapping), or D (different amino acids, not overlapping).

**Produced by:** `6-desglose.pl`

**Input:** `polyx_pair.txt` (this folder), `aa_seqs_OrthoMCL-CURRENT.fasta` (`b-Data_processed`).

**External tool:** MAFFT (pairwise alignment).

**Format:** Tab-separated, no header.

| Column | Description |
|---|---|
| 1 | Protein ID 1 |
| 2 | PolyX coordinates in protein 1 (`start-end`) |
| 3 | Amino acid in polyX of protein 1 |
| 4 | Protein ID 2 |
| 5 | PolyX coordinates in protein 2 (`start-end`) |
| 6 | Amino acid in polyX of protein 2 |
| 7 | Mean polyX length across both proteins |
| 8 | Overlap fraction (proportion of the shorter polyX aligned to the other) |
| 9 | Category (A / B / C / D) |

---

### `overlaps.txt`

**Description:** Classification of polyX pairs from ortholog pairs where at least one protein carries more than one polyX (sourced from `m1polyx.txt` and `mmpolyx.txt`). All combinations of polyX between each pair of orthologs are tested for alignment overlap. Only overlapping combinations are reported (categories A and C; non-overlapping pairs are discarded, as category B/D classification is not meaningful when one protein has multiple polyX).

**Produced by:** `7-overlap_multiple.pl`

**Input:** `m1polyx.txt` and `mmpolyx.txt` (this folder, combined into `multiple.txt` prior to running), `polyx.txt` (this folder), `aa_seqs_OrthoMCL-CURRENT.fasta` (`b-Data_processed`).

**External tool:** MAFFT (pairwise alignment).

**Format:** Tab-separated, no header.

| Column | Description |
|---|---|
| 1 | Protein ID 1 |
| 2 | PolyX coordinates in protein 1 (`start-end`) |
| 3 | Amino acid in polyX of protein 1 |
| 4 | Protein ID 2 |
| 5 | PolyX coordinates in protein 2 (`start-end`) |
| 6 | Amino acid in polyX of protein 2 |
| 7 | Mean polyX length |
| 8 | Overlap fraction |
| 9 | Category (A / C) |

**Number of entries:** 60,039.

---

## Stage 4 — Taxonomic and amino acid annotation

### `desglose_onetoone2.txt`

**Description:** Cleaned version of `desglose_onetoone.txt`, with protein IDs standardised (legacy `-old` suffixes removed) and restricted to overlapping pairs (categories A and C) used in all downstream analyses. Serves as direct input to `8-add_taxcat_aacat.pl`.

**Input:** `desglose_onetoone.txt`.

**Format:** Same nine-column tab-separated format as `desglose_onetoone.txt`.

---

### `desglose_onetoone_taxcat.txt`

**Description:** Fully annotated one-to-one pair classification table. Extends `desglose_onetoone2.txt` with two additional columns: the canonical amino acid pair code and the taxonomic category pair of the two proteins, both sorted by a predefined ranking so that each combination appears consistently.

**Produced by:** `8-add_taxcat_aacat.pl`

**Input:** `desglose_onetoone2.txt` (this folder), `tax_info.csv` (`a-Data_downloaded`).

**Format:** Tab-separated, with header. Same columns as `desglose_onetoone2.txt` plus:

| Additional column | Description |
|---|---|
| `aapair` | Canonical two-letter amino acid pair code (e.g., `AE`, `SS`) |
| `taxcat` | Taxonomic category pair (e.g., `META-FUNG`), sorted by taxonomic rank |

---

### `overlaps2.txt`

**Description:** Fully annotated multiple-polyX overlap classification table. Equivalent to `desglose_onetoone_taxcat.txt` but for pairs derived from `overlaps.txt`. Protein IDs are standardised (legacy `-old` suffixes removed) and two annotation columns are added.

**Input:** `overlaps.txt` (this folder), `tax_info.csv` (`a-Data_downloaded`).

**Format:** Tab-separated, with header.

| Column | Description |
|---|---|
| `id1` | Protein ID 1 |
| `polyx1` | PolyX coordinates in protein 1 (`start-end`) |
| `aa1` | Amino acid in polyX of protein 1 |
| `id2` | Protein ID 2 |
| `polyx2` | PolyX coordinates in protein 2 (`start-end`) |
| `aa2` | Amino acid in polyX of protein 2 |
| `mean_le` | Mean polyX length |
| `overlap` | Overlap fraction |
| `type` | Category (A / C) |
| `aapair` | Canonical amino acid pair code |
| `taxcat` | Taxonomic category pair |

---

### `desglose_multiple_taxcat.txt`

**Description:** Equivalent of `desglose_onetoone_taxcat.txt` for multiple-polyX pairs. Produced from `overlaps2.txt` by applying the same taxonomic and amino acid annotation logic. Together with `desglose_onetoone_taxcat.txt`, this file is used as combined input for all figure-generation R scripts.

**Input:** `overlaps2.txt` (this folder), `tax_info.csv` (`a-Data_downloaded`).

**Format:** Same eleven-column tab-separated format as `overlaps2.txt`.

---

## Stage 5 — Per amino acid pair files

The following files each collect all overlapping polyX pairs (categories A and C, from both `desglose_onetoone_taxcat.txt` and `desglose_multiple_taxcat.txt`) for a specific amino acid pair combination. They share a common eleven-column format (same as `overlaps2.txt`) and are used as inputs by `rel_pos.pl`, `polyAcontext.pl`, and the enrichment pipeline.

| File | Amino acid pair | Number of pairs |
|---|---|---|
| `polyAA.txt` | Alanine–Alanine (A–A) | 14,163 |
| `polyAS.txt` | Alanine–Serine (A–S) | 354 |
| `polySS.txt` | Serine–Serine (S–S) | 14,551 |
| `polyDD.txt` | Aspartate–Aspartate (D–D) | 2,614 |
| `polyDE.txt` | Aspartate–Glutamate (D–E) | 1,806 |
| `polyEE.txt` | Glutamate–Glutamate (E–E) | 15,916 |

**Format (all files):** Tab-separated, no header. Columns as in `overlaps2.txt` (id1, polyx1, aa1, id2, polyx2, aa2, mean_le, overlap, type, aapair, taxcat).

---

## Stage 6 — Relative position

### `relative.tsv`

**Description:** Relative position of every polyX region along the N→C axis of its host protein (0 = N-terminus, 1 = C-terminus), computed from the midpoint of each polyX run divided by the protein length. Both members of each pair are recorded. Covers all eight amino acid pair categories studied.

**Produced by:** `11-rel_pos.pl`

**Input:** Per-type pair files (`polyAA.txt`, `polyAS.txt`, `polySS.txt`, `polyDD.txt`, `polyDE.txt`, `polyEE.txt`; this folder), `aa_seqs_OrthoMCL-CURRENT.fasta` (`b-Data_processed`).

**Format:** Tab-separated, with header.

| Column | Description |
|---|---|
| `id` | Protein identifier |
| `relative` | Relative polyX midpoint position (0–1) |
| `polyx` | PolyX category label (e.g., `A-AA`, `E-EE`) |
| `position` | Absolute polyX coordinates (`start-end`) |

**Number of entries:** 98,820.

**PolyX category counts:** A-AA: 31,834 · S-SS: 29,104 · A-AS: 355 · S-AS: 355 · E-EE: 31,834 · D-DD: 5,230 · D-DE: 1,807 · E-DE: 1,807.

---

### `relative2.tsv`

**Description:** Filtered subset of `relative.tsv` retaining only one entry per unique protein, used for disorder analysis and for panel B of Figure 4. Deduplication ensures that proteins appearing in multiple pairs are not overrepresented in downstream statistical comparisons.

**Input:** `relative.tsv` (this folder).

**Format:** Identical four-column tab-separated format as `relative.tsv`.

**Number of entries:** 17,156.

**PolyX category counts:** E-EE: 5,248 · S-SS: 4,467 · A-AA: 4,044 · D-DD: 1,272 · E-DE: 946 · D-DE: 674 · A-AS: 261 · S-AS: 244.

---

## Stage 7 — Disorder analysis

### `polyx_disorder.txt`

**Description:** Intrinsic disorder content for each protein in `relative2.tsv`, estimated with IUPred3. For each protein, the number and fraction of residues predicted to be disordered (IUPred3 score > 0.5) are recorded alongside the polyX type and coordinates.

**Produced by:** `13-calculate_disorder.pl`

**Input:** `relative2.tsv`, `./datasets/all.fasta` (inside `enrichment.zip`).

**External tool:** IUPred3 (`iupred3.py`, short mode).

**Format:** Tab-separated, with header.

| Column | Description |
|---|---|
| `polyx` | PolyX category label (e.g., `D-DD`) |
| `id` | Protein identifier |
| `position` | PolyX coordinates (`start-end`) |
| `AAinDisorder` | Number of residues in the protein with IUPred3 score > 0.5 |
| `protlength` | Total protein length (residues) |
| `ratioINdisorder` | Fraction of residues predicted to be disordered (`AAinDisorder / protlength`) |

**Number of entries:** 17,156.

---

## Stage 8 — Sequence context and PPI model

### `A-AA.csv`

**Description:** Amino acid context (±10 flanking residues) around each polyA region in A–A pairs, extracted from the host protein sequence. Used as the prediction dataset for the Random Forest PPI model trained in `predictions.r`. The `ppi` column is set to 0 (no PPI annotation available for these project sequences).

**Produced by:** `12-polyAcontext.pl` (run with `$a = "A"`, `$aa = "A-AA"`)

**Input:** `polyAA.txt` (this folder), `filtered_seqs.fasta` (`b-Data_processed`).

**Format:** Comma-separated, with header. 24 columns: `ac`, `A-AA` (coordinates), `length`, `ppi`, followed by 20 context positions (`M10`–`M1` upstream, `P1`–`P10` downstream).

**Number of entries:** 28,328.

---

### `A-AS.csv`

**Description:** Amino acid context (±10 flanking residues) around each polyA region in A–S pairs. Same format and purpose as `A-AA.csv`. The `ppi` column is set to 0.

**Produced by:** `12-polyAcontext.pl` (run with `$a = "A"`, `$aa = "A-AS"`)

**Input:** `polyAS.txt` (this folder), `filtered_seqs.fasta` (`b-Data_processed`).

**Format:** Comma-separated, with header. Same 24-column structure as `A-AA.csv`.

**Number of entries:** 355.

---

### `polyA_model.csv`

**Description:** Labelled training dataset for the Random Forest PPI classifier implemented in `17-predictions.r`. Each row is a human polyA region with its flanking context and a binary PPI label indicating whether the protein is known to participate in protein–protein interactions. This dataset was assembled from UniProt/SwissProt human protein annotations and existing PPI databases prior to the project's main analysis.

**Format:** Comma-separated, with header. 24 columns: `ac`, `polyA` (coordinates), `length`, `ppi` (0 = no PPI, 1 = PPI), followed by 20 context positions (`M10`–`M1`, `P1`–`P10`).

**Number of entries:** 7,263 (745 PPI-positive, 6,518 PPI-negative).

**Note:** Protein IDs carry the UniProt/SwissProt `sp|` prefix (e.g., `sp|O14512|SOCS7_HUMAN`), distinguishing them from the project's OrthoMCL identifiers used in all other files.

---

## Stage 9 — GO enrichment (`enrichment.zip`)

The `enrichment.zip` archive contains all inputs and outputs for Gene Ontology enrichment analysis. It is organised into five subfolders:

---

### `enrichment.zip/annot_bg/`

**Description:** eggNOG-mapper annotation output for the background protein set, providing GO term associations used as the population background in `enrichment.py`.

**Key file:** `out.emapper.decorated.gff` — eggNOG-mapper GFF file with GO terms in the `em_GOs` attribute field of each protein entry.

**Produced by:** eggNOG-mapper, run on the background FASTA sequences (`enrichment.zip/bg/unique_id.fasta`).

---

### `enrichment.zip/bg/`

**Description:** Background sequence set for GO enrichment. Contains the protein IDs and FASTA sequences used as the population in the `goatools` enrichment analysis.

**Key files:**
- `unique_id.txt` — plain-text list of background protein IDs (one per line). Subset of the broader `unique_id.txt` in `b-Data_processed`, restricted to proteins relevant to the enrichment comparison groups. **Produced by:** `10-get_bg.pl`.
- `unique_id.fasta` — FASTA file of background protein sequences corresponding to `unique_id.txt`. **Produced by:** `10-get_bg.pl` from `filtered_seqs.fasta` (`b-Data_processed`).

---

### `enrichment.zip/datasets/`

**Description:** Per-group protein sequence and ID files used as the study sets in GO enrichment analysis.

**Key files:**
- `D-DD.fasta`, `D-DE.fasta`, `E-DE.fasta`, `E-EE.fasta` — FASTA files of proteins from each D/E amino acid pair combination, produced by `get_datasets.pl` from `desglose_onetoone_taxcat.txt` and `desglose_multiple_taxcat.txt`.
- `A-AA.fasta`, `A-AS.fasta`, `S-AS.fasta`, `S-SS.fasta` — equivalent files for the A/S group.
- `D-DD_id.txt`, `D-DE_id.txt`, `E-DE_id.txt`, `E-EE_id.txt` — plain-text gene ID lists (one ID per line) derived from the corresponding FASTA files; direct input to `enrichment.py`.
- `A-AA_id.txt`, `A-AS_id.txt`, `S-AS_id.txt`, `S-SS_id.txt` — equivalent ID lists for the A/S group.

**Produced by:** `9-get_datasets.pl` (FASTA files).

---

### `enrichment.zip/enrichment_DE/`

**Description:** GO enrichment analysis results for the D/E amino acid group (subsets D-DD, D-DE, E-DE, E-EE), produced by `enrichment.py`.

**Key files:**
- `combined_GO_enrichment_union_top5.tsv` — enrichment matrix: GO terms (rows) × subsets (columns), values are BH-corrected FDR p-values. GO terms are the union of the top 5 most enriched terms per subset.
- `GO_enrichment_union_top5_summary.tsv` — summary table with GO term ID, namespace (BP/MF/CC), and term name.
- `GO_enrichment_union_top5_grouped.png` / `.svg` — heatmap of −log10(FDR) values, grouped by GO namespace.

**Produced by:** `14-enrichment.py`

---

### `enrichment.zip/enrichment_AS/`

**Description:** GO enrichment analysis results for the A/S amino acid group (subsets A-AA, A-AS, S-AS, S-SS). Same structure and file types as `enrichment_DE/`.

**Produced by:** `14-enrichment.py` (run separately for the A/S group)

---

## Processing overview

```
b-Data_processed/filtered_pairs.txt + polyx.txt
        │
        └─► split_pairs.pl ──► 00polyx.txt, 01polyx.txt, 11polyx.txt,
                                m0polyx.txt, m1polyx.txt, mmpolyx.txt
                                        │
                    ┌───────────────────┴──────────────────────┐
                    │                                           │
              get_info.pl                          overlap_multiple.pl
                    │                                           │
             polyx_pair.txt                           overlaps.txt
                    │                                           │
              desglose.pl                          [annotation + cleaning]
                    │                                           │
          desglose_onetoone.txt                       overlaps2.txt
                    │                                           │
          [clean → desglose_onetoone2.txt]                      │
                    │                                           │
          add_taxcat_aacat.pl                   desglose_multiple_taxcat.txt
                    │                                           │
        desglose_onetoone_taxcat.txt                            │
                    └───────────────────┬───────────────────────┘
                                        │
                            [subset by amino acid pair]
                                        │
               polyAA.txt, polyAS.txt, polySS.txt,
               polyDD.txt, polyDE.txt, polyEE.txt
                                        │
                    ┌───────────────────┼──────────────────────┐
                    │                   │                       │
               rel_pos.pl       polyAcontext.pl         get_datasets.pl
                    │                   │                       │
            relative.tsv       A-AA.csv, A-AS.csv       datasets/ (FASTA + IDs)
                    │                   │                       │
            [deduplicate]       predictions.r            enrichment.py
                    │                                           │
            relative2.tsv                              enrichment_DE/, enrichment_AS/
                    │
        calculate_disorder.pl
                    │
          polyx_disorder.txt
```
