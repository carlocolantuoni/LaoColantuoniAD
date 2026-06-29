# NPTX2-Centered Cognitive Resilience Mechanisms in the Context of AD Pathology

Analysis and figure-generation code accompanying:

> Lao et al. (2026). *NPTX2-Centered Cognitive Resilience Mechanisms in the Context
> of AD Pathology.* DOI: [10.1101/2025.10.17.683150](https://doi.org/10.1101/2025.10.17.683150)

The pipeline characterizes NPTX2 co-expression/co-abundance trajectories across the
Alzheimer's disease continuum — comparing pathology-defined control subgroups
(CN-Lo, CN-Hi) against MCI and AD, in both RNA and protein (PRM-MS) modalities —
and produces every results table and figure panel reported in the paper.

> **License:** *(add a license, e.g. MIT or CC-BY-4.0; see "License" below)*

## Repository contents

| File | Purpose |
| --- | --- |
| `01_analysis.R` | All computation, no plotting. Reads the curated input data and writes every results table. |
| `02_visualization.R` | All plotting, no computation. Reads the tables from `01` and renders the figures (PNG + SVG). |

Run `01_analysis.R` first; it produces the tables that `02_visualization.R` consumes.

## Requirements

Developed and run under R (>= 4.1 recommended).

CRAN packages (required):

- `dplyr`, `readr`, `tidyr`, `tibble`
- `ggplot2`, `ggpubr`, `cowplot`, `svglite`

Bioconductor / optional packages — used by individual analysis steps, which are
**skipped gracefully with a message** if the package is absent, so the rest of the
pipeline still runs:

- `limma`, `matrixStats` (differential expression, Supp 4a/4b)
- `clusterProfiler`, `org.Hs.eg.db` (over-representation analysis, Fig 2f/3g)
- `fgsea` (gene-set enrichment, Supp 2b/2c)
- `biomaRt` (protein-coding gene filter; only queried if a local symbol list is not provided — see below)

Install everything with:

```r
install.packages(c("dplyr", "readr", "tidyr", "tibble",
                   "ggplot2", "ggpubr", "cowplot", "svglite"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "matrixStats", "clusterProfiler",
                       "org.Hs.eg.db", "fgsea", "biomaRt"))
```

## Input data

`01_analysis.R` reads three files from an input folder (default `exported_data_fixed/`):

| File | Description |
| --- | --- |
| `metadata_master.csv` | Per-sample metadata (group, diagnosis, sex, age, PMI, RIN, Braak/CERAD, PC1/PC2, PRM membership). |
| `protein_abundance.csv` | Per-sample protein measurements: `MS_*` (PRM-MS) and `WB_*` columns. |
| `expression_rna.csv.gz` | RNA expression matrix (genes x samples; first column is the gene symbol). |

Optional, in the same folder:

- `protein_coding_symbols.csv` — a one-column list of protein-coding gene symbols.
  If present it is used for the protein-coding filter (fully reproducible, no network);
  if absent, the code falls back to a `biomaRt` query, and if that is unavailable the
  filter is skipped.

The input data are not included in this repository. *(Add a sentence here pointing to
where they can be obtained — e.g. a data repository, accession number, or
"available from the authors on request.")*

### Gene-set file (for the fGSEA step)

The fGSEA step (Supp 2b/2c) needs a GO:CC `.gmt` file (e.g. an MSigDB C5 GO:CC
collection such as `c5.go.cc.v2026.1.Hs.symbols.gmt`). Place it in a `gene_sets/`
folder at the repository root, or point to it with an environment variable:

```bash
export GMT_DIR=/path/to/your/gmt/folder
```

If no `.gmt` file is found, this single step is skipped (the PC1/PC2 loading-rank
vectors are still written); the rest of the pipeline is unaffected.

## Running the pipeline

From an R session with the working directory set to the repository root:

```r
source("01_analysis.R")        # writes result tables
source("02_visualization.R")   # reads those tables, writes figures
```

Or from the shell:

```bash
Rscript 01_analysis.R
Rscript 02_visualization.R
```

### Selecting the input/output set

`01_analysis.R` has a single switch near the top, `INPUT_FOLDER`. The output folder
is derived from it automatically so separate runs never overwrite each other:

- `INPUT_FOLDER <- "exported_data"`        -> results in `analysis_output/`
- `INPUT_FOLDER <- "exported_data_fixed"`  -> results in `analysis_output_fixed/`

`02_visualization.R` has matching `DATA_FOLDER` / `RESULTS_FOLDER` / `FIG_FOLDER`
settings at its top; keep these consistent with the analysis run you want to plot.

## Outputs

`01_analysis.R` -> result tables (`.csv`) in `analysis_output_fixed/`, one set per
figure panel (boxplot values and comparison statistics, trajectory metrics,
correlation distributions, the master gene classification, double differential
expression, sex diagnostics, Braak/CERAD composition, and the supplementary tables).

`02_visualization.R` -> figures in `final_figures_fixed/`, written as paired
**PNG** (`png/`) and **SVG** (`svg/`) for every panel, covering main Figures 1–4
and the supplementary figures.

## Notes

- **Group labels.** Internally the pathology subgroups are relabeled
  `LoPath -> CN-Lo` and `HiPath -> CN-Hi` (sample membership unchanged); the CN-Hi
  prefix carries through to the trajectory-pattern names.
- **Reproducibility.** Stochastic steps (fGSEA) set a fixed seed. Enrichment results
  can still shift slightly with the installed `org.Hs.eg.db` / `GO.db` version.

## License

*(Choose and add one — e.g. drop a `LICENSE` file in the repo. MIT is common for
code; CC-BY-4.0 is common when code accompanies a publication.)*

## Contact

*(Corresponding author name and email, or a link to open an issue.)*
