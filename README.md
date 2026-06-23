# b-Data_processed — Processed Data

This folder contains data files derived from the downloaded data in `a-Data_downloaded` through filtering and extraction steps. 

---

## Files

### 1. `filtered_pairs.zip`

**Description:** Compressed archive containing `filtered_pairs.txt`, a subset of the full ortholog pair table retaining only pairs whose normalized similarity score exceeds the median score across all pairs. This filtering step removes low-confidence ortholog relationships before all downstream polyX analyses.

**Produced by:** `2-filter_pairs.pl` (see `d-Software/README.md`)

**Input:** `orthologs.txt` (from `orthologs.zip`, folder `a-Data_downloaded`)

**Format (`filtered_pairs.txt`, inside the archive):** Tab-separated, inheriting the column structure of the source `orthologs.txt`:

| Column | Description |
|---|---|
| 1 | Protein ID 1 |
| 2 | Protein ID 2 |
| 3 | Normalized similarity score (> median) |

**Modifications relative to source:** Rows with a normalized similarity score ≤ median were removed. No other changes to content or format.

---

### 2. `unique_id.txt`

**Description:** Plain-text list of unique protein IDs present in `filtered_pairs.txt`. Each ID appears exactly once, regardless of how many ortholog pairs it participates in. This file serves as the ID filter for sequence extraction and as the background gene list for GO enrichment analysis.

**Produced by:** Extracted from `2-filtered_pairs.txt` by collecting all values in columns 1 and 2 and retaining unique entries.

**Input:** `filtered_pairs.txt` (from `filtered_pairs.zip`, this folder)

**Format:** Plain text, one protein ID per line, no header.

**Modifications relative to source:** Derived from `filtered_pairs.txt` by deduplication; no sequence or annotation data included.

---
