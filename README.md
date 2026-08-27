# WntQuant

Directional quantification of Wnt signaling-pathway activity from bulk RNA-seq
data.

## What is WntQuant?

WntQuant derives and refines gene signatures of Wnt pathway activity from
multiple expression datasets. Given samples split into high-Wnt (H) and
low-Wnt (L) groups, it returns two kinds of genes:

- **Activation genes**: higher expression in the high-Wnt group.
- **Inhibition genes**: higher expression in the low-Wnt group.

It then removes redundant gene sets and returns compact, high-confidence
signatures.

## How it works

WntQuant has three steps, each implemented as one function:

1. `Get_Wnt_denovo_genesets` - runs differential analysis (limma, t-test, or
   Wilcoxon test) between the H and L groups and ranks genes into activation
   and inhibition signatures.
2. `Wnt_purification_system` - computes a "Score2" confidence metric to clean
   de novo signatures, or validates external Wnt gene sets using fGSEA.
3. `Merge_Wnt_genesets` - drops gene sets that are too small and merges highly
   similar ones using Jaccard similarity.

## Installation

Run the following in a fresh R session:

```r
# CRAN dependencies
install.packages("BiocManager")
install.packages(c("doBy", "dplyr", "igraph", "magrittr", "metap", "reshape2"))

# Bioconductor dependencies
BiocManager::install(c("limma", "fgsea", "impute", "survcomp"))

# WntQuant itself
install.packages("remotes")
remotes::install_github("FangZY-Lab/WntQuant")
```

Then load it:

```r
library(WntQuant)
```

## Data format

WntQuant reads data from the global environment via `get()`. For each dataset
you need:

- An expression data frame named `<Dataset>` with genes as rows and samples as
  columns.
- A group data frame named `<Dataset>_G` with a `Tag` column (sample IDs
  matching the expression column names) and a `group` column (group labels).
- A `group_HL` data frame with columns `Accession`, `Group1`, `Group2`,
  `Group1_Status`, and `Group2_Status`, where each status is `"H"` or `"L"`.

`file_paths` is the directory where result files are written when
`export_file = TRUE`.

## Quick start example

This example creates two simulated datasets and runs the full pipeline:

```r
library(WntQuant)
set.seed(123)

make_dataset <- function(n_genes = 100, n_h = 30, n_l = 30, seed = 1) {
  set.seed(seed)
  n_samples <- n_h + n_l
  samples <- c(paste0("S", seq_len(n_h)), paste0("S", n_h + seq_len(n_l)))

  mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)
  # Activation genes: higher in the H group.
  mat[1:10, seq_len(n_h)] <- mat[1:10, seq_len(n_h)] + 2
  # Inhibition genes: higher in the L group.
  mat[11:20, n_h + seq_len(n_l)] <- mat[11:20, n_h + seq_len(n_l)] + 2

  expr <- as.data.frame(mat)
  rownames(expr) <- paste0("GENE", seq_len(n_genes))
  colnames(expr) <- samples

  group <- data.frame(
    Tag = samples,
    group = c(rep("GroupA", n_h), rep("GroupB", n_l)),
    stringsAsFactors = FALSE
  )
  list(expr = expr, group = group)
}

dA <- make_dataset(seed = 11)
dB <- make_dataset(seed = 22)

SimDataA <- dA$expr
SimDataA_G <- dA$group
SimDataB <- dB$expr
SimDataB_G <- dB$group

group_HL <- data.frame(
  Accession = c("SimDataA", "SimDataB"),
  Group1 = c("GroupA", "GroupA"),
  Group2 = c("GroupB", "GroupB"),
  Group1_Status = c("H", "H"),
  Group2_Status = c("L", "L"),
  stringsAsFactors = FALSE
)

out_dir <- tempdir()

# 1. De novo activation / inhibition signatures.
denovo <- Get_Wnt_denovo_genesets(
  file_paths = out_dir,
  expression_accession_vector = c("SimDataA", "SimDataB"),
  group_HL = group_HL,
  gene_difference_method = "limma",
  alternative = "two.sided",
  threshold_or_rank = "rank",
  top_genes = 30
)

# 2. Purify the signatures with the Score2 metric.
refined <- Wnt_purification_system(
  file_paths = out_dir,
  expression_accession_vector = c("SimDataA", "SimDataB"),
  group_HL = group_HL,
  activation_geneset = denovo$activation,
  inhibition_geneset = denovo$inhibition,
  purpose = "cleaned",
  threshold_type = "quantile",
  quantile_threshold = 0.9
)

# 3. Merge redundant signatures.
merged <- Merge_Wnt_genesets(
  file_paths = out_dir,
  activation_geneset = refined$update_activation_geneset,
  inhibition_geneset = refined$update_inhibition_geneset,
  delete_GN = TRUE,
  min_GN = 3
)

str(merged$joint_activation_geneset)
str(merged$joint_inhibition_geneset)
```

If it works, the final activation signature contains the planted activation
genes (`GENE1`-`GENE10`) and the inhibition signature contains the planted
inhibition genes (`GENE11`-`GENE20`).

## Functions

| Function | Purpose |
| --- | --- |
| `Get_Wnt_denovo_genesets()` | Differential analysis plus ranking to derive activation / inhibition signatures |
| `Wnt_purification_system()` | Score2-based cleaning or fGSEA validation |
| `Merge_Wnt_genesets()` | Size filtering plus Jaccard-based merging |

## Notes

1. Data must be in the global environment and named exactly `<Dataset>` and
   `<Dataset>_G`.
2. Use at least two datasets, as in the example.
3. `file_paths` is only used for writing result files when `export_file =
   TRUE`.

## Author

Dingkang Zhao (赵定康)

Email: <dingkang.25@intl.zju.edu.cn>
