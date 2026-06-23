# d-Software — Code Repository

This folder contains all software generated during the project for processing and analyzing polyX regions in orthologous protein sequences across eukaryotic and prokaryotic taxa. Scripts are written in Perl, Python, and R.

Each script is described below with its purpose, inputs, and outputs.

---

## Table of Contents

1. [1-polyx2_standalone.pl]
2. [2-filter_pairs.pl]
3. [3-filter_fasta.pl]
4. [4-split_pairs.pl]
5. [5-get_info.pl]
6. [6-desglose.pl]
7. [7-overlap_multiple.pl]
8. [8-add_taxcat_aacat.pl]
9. [9-get_datasets.pl]
10. [10-get_bg.pl]
11. [11-rel_pos.pl]
12. [12-polyAcontext.pl]
13. [13-calculate_disorder.pl]
14. [14-enrichment.py]
15. [15-plots.R]
16. [16-relative.R]
17. [17-predictions.r]
18. [18-fig4.R]
19. [19-fig5.R]

---

## Perl Scripts

### 1. `polyx2_standalone.pl`

**Description:** Core polyX scanner. Scans protein sequences from a FASTA file and identifies polyX regions using a sliding-window threshold (default: ≥8 identical residues in a window of 10). Overlapping windows are merged into a single run. The beginning and end of each run are trimmed of non-specific residues.

**Usage:**
```
perl 1-polyx2_standalone.pl <INPUT.fasta> [MINIMUM] [WINDOW]
```
- `MINIMUM`: minimum number of the same amino acid in the window (default: 8)
- `WINDOW`: window length (default: 10)

**Input:**
- A protein sequence file in FASTA format.

**Output:**
- `polyx.txt` — tab-separated file with columns: `ID`, `Start`, `End`, `Aa`, `polyX` (the sequence of the run).

**Dependencies:** BioPerl (`Bio::SeqIO`)

---

### 2. `filter_pairs.pl`

**Description:** Filters ortholog pairs by normalized similarity score. Computes the median score across all pairs and retains only those above it.

**Input:**
- `orthologs.txt` — tab-separated file of ortholog pairs; the third column is the normalized similarity score.

**Output:**
- `filtered_pairs.txt` — subset of the input containing only above-median pairs.

---

### 3. `filter_fasta.pl`

**Description:** Extracts a subset of sequences from a FASTA file given a list of protein IDs.

**Input:**
- `unique_id.txt` — plain-text file with one protein ID per line.
- `aa_seqs_OrthoMCL-CURRENT.fasta` — full protein sequence database in FASTA format.

**Output:**
- `filtered_seqs.fasta` — FASTA file containing only the sequences whose IDs appear in `unique_id.txt`.

---

### 4. `split_pairs.pl`

**Description:** Classifies ortholog pairs based on how many polyX regions each member carries (0, 1, or >1), and writes each class to a separate file.

**Input:**
- `polyx.txt` — polyX table produced by `1-polyx2_standalone.pl`.
- `filtered_pairs.txt` — filtered ortholog pairs from `2-filter_pairs.pl`.

**Output:**
- `00polyx.txt` — both proteins have 0 polyX.
- `01polyx.txt` — one protein has 0, the other has 1 polyX.
- `11polyx.txt` — both proteins have exactly 1 polyX.
- `m1polyx.txt` — one protein has >1 polyX, the other has exactly 1.
- `mmpolyx.txt` — both proteins have >1 polyX.
- `m0polyx.txt` — one protein has >1 polyX, the other has 0.

---

### 5. `get_info.pl`

**Description:** Annotates ortholog pairs with polyX coordinates and amino acid type for each member by looking up data from the polyX table.

**Input:**
- `11polyx.txt` — ortholog pairs where both proteins have exactly 1 polyX (from `4-split_pairs.pl`).
- `polyx.txt` — polyX coordinates and amino acid type per protein.

**Output:**
- `polyx_pair.txt` — tab-separated file with columns: `ID1`, `Aa1`, `Start1`, `End1`, `ID2`, `Aa2`, `Start2`, `End2`.

---

### 6. `desglose.pl`

**Description:** Classifies one-to-one ortholog polyX pairs (both members have exactly 1 polyX) into four categories based on pairwise sequence alignment (via MAFFT) and the overlap of polyX regions:
- **A** (`same_overlap`): same amino acid type, polyX regions overlap.
- **B** (`same_notaligned`): same amino acid type, regions do not overlap.
- **C** (`diff_aligned`): different amino acid types, regions overlap.
- **D** (`diff_notaligned`): different amino acid types, regions do not overlap.

Overlap is defined as ≥50% of the shorter polyX aligned to the other.

**Input:**
- `polyx_pair.txt` — produced by `5-get_info.pl`.
- `aa_seqs_OrthoMCL-CURRENT.fasta` — full protein sequence database.

**Output:**
- `desglose_onetoone.txt` — tab-separated file with columns: `ID1`, `Range1`, `Aa1`, `ID2`, `Range2`, `Aa2`, `Mean_length`, `Overlap_fraction`, `Category`.

**External tools:** MAFFT (must be installed and in PATH).

**Dependencies:** BioPerl (`Bio::SeqIO`)

---

### 7. `overlap_multiple.pl`

**Description:** Processes pairs where at least one protein has more than one polyX region. All combinations of polyX between the two proteins are tested for alignment overlap (via MAFFT), and overlapping pairs are classified as:
- **A**: same amino acid type and ≥50% overlap.
- **C**: different amino acid types and ≥50% overlap.

**Input:**
- `multiple.txt` — ortholog pairs where at least one protein has >1 polyX.
- `aa_seqs_OrthoMCL-CURRENT.fasta` — full protein sequence database.
- `polyx.txt` — polyX coordinates per protein.

**Output:**
- `overlaps.txt` — tab-separated file with the same column structure as `desglose_onetoone.txt` (ID1, Range1, Aa1, ID2, Range2, Aa2, Mean_length, Overlap_fraction, Category).

**External tools:** MAFFT (must be installed and in PATH).

**Dependencies:** BioPerl (`Bio::SeqIO`), `File::Temp`

---

### 8. `add_taxcat_aacat.pl`

**Description:** Annotates a desglose file with two additional columns: (1) the canonical amino acid pair (`aapair`), sorted alphabetically by a predefined amino acid ranking; and (2) the taxonomic category pair (`taxcat`), sorted by a predefined taxonomic ranking using species abbreviations from a reference table.

**Input:**
- `desglose_onetoone2.txt` — annotated one-to-one ortholog pairs (hardcoded).
- `tax_info.csv` — comma-separated table mapping species IDs to taxonomy category abbreviations.

**Output:**
- `desglose_onetoone_taxcat.txt` — input file with two appended columns: `aapair` and `taxcat`.

**Dependencies:** BioPerl (`Bio::SeqIO`)

---

### 9. `get_datasets.pl`

**Description:** Extracts protein sequences into a FASTA file for a specific amino acid pair (e.g., `EE`) from the overlap classifications produced by both `6-desglose.pl` (one-to-one) and `7-overlap_multiple.pl` (multiple). Only proteins involved in overlapping pairs (categories A or C) matching the target amino acid pair are included, with deduplication.

**Input:**
- `filtered_seqs.fasta` — filtered protein sequence database.
- `desglose_onetoone.txt` — one-to-one classification output.
- `overlaps.txt` — multiple polyX classification output.

**Output:**
- `<aa1>-<pair>.fasta` — FASTA file of sequences for the specified amino acid pair (e.g., `E-EE.fasta`).

**Note:** The target amino acids (`aa1`, `aa2`, `pair`) must be set by editing the script variables directly.

---

### 10. `get_bg.pl`

**Description:** Extracts a background set of protein sequences from the filtered FASTA database, given a list of protein IDs. Used to prepare the background for GO enrichment analysis.

**Input:**
- `unique_id.txt` — plain-text file with one protein ID per line.
- `filtered_seqs.fasta` — filtered protein sequence database.

**Output:**
- `unique_id.fasta` — FASTA file of background sequences.

---

### 11. `rel_pos.pl`

**Description:** Computes the relative position of each polyX region within its host protein (from 0 = N-terminus to 1 = C-terminus), based on the midpoint of the polyX coordinates and the total protein length. Both proteins in each pair are processed.

**Input:**
- `polyEE.txt` — polyX pair file for a specific amino acid type (filename must be edited for other types).
- `aa_seqs_OrthoMCL-CURRENT.fasta` — full protein sequence database (used to retrieve sequence lengths).

**Output:**
- `relative.tsv` — tab-separated file with columns: `id`, `relative_position`, `polyx_type`, `coordinates`, `taxonomic_category`.

---

### 12. `polyAcontext.pl`

**Description:** Extracts the local amino acid context (±10 residues) flanking the start and end of each polyX region for a given amino acid type. Processes both proteins in each pair, with deduplication. Unknown or ambiguous residues (`X`) are replaced by `-`.

**Input:**
- `polyAA.txt` — pair file for the target amino acid type (hardcoded; edit `$a` and `$aa` variables for other types).
- `filtered_seqs.fasta` — filtered protein sequence database.

**Output:**
- `<aa>-<pair>.csv` (e.g., `A-AA.csv`) — comma-separated file with columns: `ac`, coordinates, `length`, `ppi`, and 20 context positions (`M10` to `M1` upstream, `P1` to `P10` downstream).

---

### 13. `calculate_disorder.pl`

**Description:** Estimates the intrinsic disorder content of proteins harboring polyX regions using IUPred3. For each protein, the fraction of residues with an IUPred3 score >0.5 (disordered) is calculated.

**Input:**
- `all.fasta` — FASTA file of all relevant protein sequences.
- `relative2.tsv` — table of relative polyX positions per protein and type.

**Output:**
- `polyx_disorder.txt` — tab-separated file with columns: `polyx` (type), `id`, `position`, `AAinDisorder` (count of disordered residues), `protlength`, `ratioINdisorder`.

**External tools:** IUPred3 (`./IUPred/iupred3.py`, must be installed locally); Python 3.

**Dependencies:** BioPerl (`Bio::SeqIO`)

---

## Python Scripts

### 14. `enrichment.py`

**Description:** Performs Gene Ontology (GO) enrichment analysis for four protein subsets (D-DD, D-DE, E-DE, E-EE) using `goatools`. Background annotations are loaded from an eggNOG-mapper GFF file. The top 5 most enriched GO terms per subset are selected, and their union is used to build a combined enrichment matrix. Results are visualized as a heatmap grouped by GO namespace (MF, CC, BP), with a namespace color bar.

**Input:**
- `./annot_bg/out.emapper.decorated.gff` — eggNOG-mapper annotated GFF file with GO term associations for background genes.
- `./datasets/D-DD_id.txt`, `D-DE_id.txt`, `E-DE_id.txt`, `E-EE_id.txt` — plain-text files with one gene ID per line, one per subset.
- `go-basic.obo` — GO ontology file in OBO format.

**Output (in `./enrichment_DE/`):**
- `combined_GO_enrichment_union_top5.tsv` — full enrichment matrix (GO terms × subsets, FDR-corrected p-values).
- `GO_enrichment_union_top5_summary.tsv` — summary table with GO ID, namespace, and term name.
- `GO_enrichment_union_top5_grouped.png` and `.svg` — heatmap figure (−log10 FDR, grouped by namespace).

**Dependencies:** `pandas`, `goatools`, `matplotlib`, `numpy`

---

## R Scripts

### 15. `plots.R`

**Description:** Generates three figures from the combined one-to-one and multiple polyX classification tables (annotated with taxonomy and amino acid pair information):
- **SupplFig1** — density and violin plots comparing mean polyX length between same-type (A) and different-type (C) overlapping pairs.
- **Fig2** — three-panel figure: (A) ECDF of polyX overlap fraction by type; (B) stacked bar chart of type proportions per taxonomic category pair; (C) stacked bar chart of amino acid composition per type.
- **Fig3** — heatmaps of normalized amino acid pair counts per taxonomic category, separately for type A (same polyX type) and type C (different polyX type).

**Input:**
- `desglose_onetoone_taxcat.txt` — annotated one-to-one pairs (from `8-add_taxcat_aacat.pl`).
- `desglose_multiple_taxcat.txt` — annotated multiple pairs (equivalent annotation for multiple-polyX pairs).

**Output:**
- `SupplFig1.pdf`
- `Fig2.pdf`
- `Fig3.pdf`

**Dependencies:** `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `ggpubr`, `FSA`, `ggnewscale`, `RColorBrewer`

---

### 16. `relative.R`

**Description:** Generates a notched boxplot of polyX relative position (N→C axis, 0–1) across eight polyX categories (A-AA, A-AS, S-AS, S-SS, D-DD, D-DE, E-DE, E-EE), grouped by amino acid family (A/S in orange, D/E in purple). Pairwise Wilcoxon tests are shown for within-family comparisons. Sample sizes are displayed below each box.

**Input:**
- `relative2.tsv` — tab-separated file with columns `id`, `relative` (position), `polyx` (category), `coordinates`, and `taxa`; produced by `11-rel_pos.pl`.

**Output:**
- `Fig4a.pdf` — boxplot figure.

**Dependencies:** `ggplot2`, `ggpubr`

---

### 17. `predictions.r`

**Description:** Trains a Random Forest classifier (via `caret`, with repeated 10-fold cross-validation and downsampling) to predict protein–protein interaction (PPI) involvement from polyX context features. The trained model is then applied to two new datasets (A-AA and A-AS) to generate PPI probability scores. Results are visualized as a boxplot (score distribution per polyX type, with Wilcoxon test) and a scatter plot (score vs. polyX length, colored by type).

**Input:**
- `./polyA_model.csv` — training dataset with context features and a binary `ppi` label (0/1); produced from the output of `12-polyAcontext.pl` with added PPI annotation.
- `A-AA.csv`, `A-AS.csv` — feature files for prediction, in the same format as the training data (from `12-polyAcontext.pl`).

**Output:**
- `score.pdf` — two-panel figure: boxplot of predicted PPI scores and scatter plot of score vs. length.

**Dependencies:** `readr`, `caret`, `pROC`, `dplyr`, `gridExtra`, `ggplot2`, `ggpubr`

---

### 18. `fig4.R`

**Description:** Generates a two-panel figure:
- **Panel A** — scatter plot with linear regression showing the relationship between the number of different proteins and the number of polyX in overlapping regions across amino acid pairs, with R² annotation.
- **Panel B** — notched boxplot of relative polyX position (N→C) across eight polyX categories, grouped by amino acid family, with pairwise Wilcoxon tests.

**Input:**
- `diffprot.txt` — table with columns `prot` (number of different proteins) and `polyx` (number of polyX in overlapping regions), one row per amino acid pair.
- `relative2.tsv` — relative polyX position table (from `11-rel_pos.pl`).

**Output:**
- `Fig4.pdf` — two-panel figure.

**Dependencies:** `ggplot2`, `ggpubr`, `ggrepel`, `gridExtra`, `reshape2`, `RColorBrewer`, `ape`

---

### 19. `fig5.R`

**Description:** Generates a two-panel GO enrichment heatmap figure from the results of `14-enrichment.py`, for two amino acid group comparisons (D-E and A-S). Each panel shows a heatmap of −log10(FDR) values for GO terms (rows) across four protein subsets (columns), with a color-coded namespace ribbon (MF/CC/BP) on the right. GO terms are grouped and sorted by namespace and mean enrichment.

**Input:**
- `./enrichment_DE/GO_enrichment_union_top5_summary.tsv` — GO term summary for the D-E group.
- `./enrichment_DE/combined_GO_enrichment_union_top5.tsv` — enrichment matrix for the D-E group.
- `./enrichment_AS/GO_enrichment_union_top5_summary.tsv` — GO term summary for the A-S group.
- `./enrichment_AS/combined_GO_enrichment_union_top5.tsv` — enrichment matrix for the A-S group.

**Output:**
- `Fig5-2.pdf` — combined two-panel figure.

**Dependencies:** `ggplot2`, `dplyr`, `readr`, `tidyr`, `patchwork`, `forcats`, `stringr`

---

## General Requirements

| Language | Version tested | Key dependencies |
|----------|---------------|-----------------|
| Perl | ≥5.26 | BioPerl (`Bio::SeqIO`), `File::Temp` |
| Python | ≥3.8 | `pandas`, `goatools`, `matplotlib`, `numpy` |
| R | ≥4.1 | `ggplot2`, `ggpubr`, `caret`, `pROC`, `patchwork`, and others listed per script |

External tools required: **MAFFT** (sequence alignment, used by `6-desglose.pl` and `7-overlap_multiple.pl`) and **IUPred3** (disorder prediction, used by `13-calculate_disorder.pl`).
