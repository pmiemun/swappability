# a-Data_downloaded — Downloaded Data

This folder contains all external data downloaded for the project prior to any processing. Data originates from public databases and ontology resources.

---

## Files

### 1. `raw_species.ods`

**Description:** Master spreadsheet cataloguing all species included in the study. Each row corresponds to one species (or strain), with its taxonomic classification, database of origin, proteome version, and sequence statistics.

**Source:** Multiple VEuPathDB-associated databases and UniProt, accessed individually per species:

**Number of species:** 865

**Columns:**

| Column | Description |
|---|---|
| `Category` | Taxonomic category abbreviation (see `taxonomic_categories.txt`) |
| `Name` | Full species/strain name |
| `TaxID` | NCBI Taxonomy ID |
| `Core/Peripheral` | Whether the species belongs to the core or peripheral set of the OrthoMCL analysis |
| `Abbreviation` | Short species code used throughout the project (see `tax_info.csv`) |
| `Resource` | Database from which the proteome was downloaded |
| `Proteome_Version` | Version or build identifier of the proteome at time of download |
| `Most_recent_version` | Whether the downloaded version was the most recent available (`most_recent` or a date string) |
| `Sequences` | Total number of sequences in the downloaded proteome |
| `Grouped_sequences` | Number of sequences included in OrthoMCL groups |
| `Groups` | Number of OrthoMCL groups the species contributes to |

**Modifications:** This file was created and maintained manually during the project to track the provenance of all downloaded proteomes.

---

### 2. `tax_info.csv`

**Description:** Reference table mapping species abbreviations to their taxonomic category codes.

**Source:** Derived from `raw_species.ods` (columns `Category` and `Abbreviation`).

**Format:** Comma-separated, no index. Two columns:

| Column | Description |
|---|---|
| `cat` | Taxonomic category abbreviation (e.g., `FUNG`, `META`) |
| `abbrev` | Species abbreviation (e.g., `bbig`, `aaeg`) |

---

### 3. `taxonomic_categories.txt`

**Description:** Plain-text legend defining the ten taxonomic category abbreviations used throughout the project, mapping each code to its full biological group name.

**Format:** Plain text, one entry per line, in the format `ABBREVIATION => Full name`.

**Contents:**

| Abbreviation | Group |
|---|---|
| ALVE | Alveolata |
| AMOE | Amoeba |
| ARCH | Archaea |
| EUGL | Euglenozoa |
| FUNG | Fungi |
| META | Metazoa |
| OBAC | Other Bacteria |
| OEUK | Other Eukaryota |
| PROT | Proteobacteria |
| VIRI | Viridiplantae |

---

### 4. `orthologs.zip`

**Description:** Compressed archive containing the ortholog pairs used as the starting point for all downstream analyses. The file `orthologs.txt` inside the archive lists pairwise orthologous relationships between proteins from different species, together with a normalized similarity score for each pair.

**Source:** OrthoMCL Database (OrthoMCL-DB), downloaded from [https://orthomcl.org](https://orthomcl.org).

**Version / release:** 6.9

**Format (orthologs.txt, inside the archive):** Tab-separated, three columns: protein ID 1, protein ID 2, normalized similarity score.

---

### 5. `go-basic.obo`

**Description:** The Gene Ontology (GO) database in OBO flat-file format, restricted to the basic subset. This file defines the GO term hierarchy (identifiers, names, namespaces, and parent–child relationships) across the three GO domains: Biological Process (BP), Molecular Function (MF), and Cellular Component (CC).

**Source:** Gene Ontology Consortium, downloaded from [https://geneontology.org/docs/download-ontology/](https://geneontology.org/docs/download-ontology/) (`go-basic.obo` — the filtered, non-redundant version recommended for most enrichment analyses).

**Date of download:** 07th June 2025

**GO release date (internal header):** 2025-06-01

---
