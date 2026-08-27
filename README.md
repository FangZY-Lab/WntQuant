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

WntQuant reads data from the calling environment via `get()`. Create the data
objects in the same environment where you call the functions (for example, at
the top level or in the same R Markdown chunk). For each dataset you need:

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

## Function reference

### `Get_Wnt_denovo_genesets()`

Differential analysis plus ranking to derive activation / inhibition
signatures.

```r
Get_Wnt_denovo_genesets(
  file_paths = out_dir,
  expression_accession_vector = c("SimDataA", "SimDataB"),
  group_HL = group_HL,
  gene_difference_method = "limma",
  alternative = "two.sided",
  p_combine_method = "fisher",
  threshold_or_rank = "rank",
  top_genes = 500,
  gene_Pfilter = 0.05,
  gene_FCfilter = 1,
  export_file = FALSE
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `file_paths` | Working directory where result files are written when `export_file = TRUE`. | required |
| `expression_accession_vector` | Character vector of dataset names. Each must exist in the global environment as an expression data frame (genes as rows, samples as columns), with a matching `<Dataset>_G` group data frame. | required |
| `group_HL` | Data frame that maps group labels to Wnt high/low status. Columns: `Accession`, `Group1`, `Group2`, `Group1_Status`, `Group2_Status`, where each status is `"H"` or `"L"`. | required |
| `gene_difference_method` | Differential test used between the H and L groups: `"limma"`, `"t_test"`, or `"wilcox_test"`. | `"limma"` |
| `alternative` | Test direction: `"two.sided"`, `"less"`, or `"greater"`. | `"two.sided"` |
| `p_combine_method` | p-value integration across sub-datasets sharing the same prefix: `"fisher"`, `"z.transform"`, `"logit"`, `"cct"`, `"sumz"`, or `"geometric_mean"`. | `"fisher"` |
| `threshold_or_rank` | Gene selection strategy: `"rank"` keeps the top-ranked genes by p-value, `"threshold"` uses p-value and fold-change cutoffs. | `"rank"` |
| `top_genes` | Number of top-ranked genes retained when `threshold_or_rank = "rank"`. | `500` |
| `gene_Pfilter` | p-value cutoff when `threshold_or_rank = "threshold"`. | `0.05` |
| `gene_FCfilter` | log2 fold-change cutoff when `threshold_or_rank = "threshold"`. | `1` |
| `export_file` | Logical. If `TRUE`, write the activation / inhibition tables to tab-delimited files in `file_paths`. | `FALSE` |

Returns a list with `activation` and `inhibition` data frames (one column per
dataset prefix, containing ranked gene IDs).

### `Wnt_purification_system()`

Computes a "Score2" confidence metric to clean de novo signatures, or validates
external Wnt gene sets with fGSEA.

```r
Wnt_purification_system(
  file_paths = out_dir,
  expression_accession_vector = c("SimDataA", "SimDataB"),
  group_HL = group_HL,
  gene_difference_method = "limma",
  alternative = "two.sided",
  p_combine_method = "fisher",
  using_FC = FALSE,
  na_ratio = 0.5,
  using_KNN = TRUE,
  statistics = "arithmetic_mean",
  purpose = "cleaned",
  threshold_type = "quantile",
  rank_threshold = 100,
  quantile_threshold = 0.99,
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = NA,
  min.sz = 1,
  max.sz = 10000,
  export_file = FALSE
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `file_paths` | Working directory where result files are written when `export_file = TRUE`. | required |
| `expression_accession_vector` | Character vector of dataset names (same requirement as above). | required |
| `group_HL` | Group-to-H/L mapping (same format as above). | required |
| `gene_difference_method` | Differential test: `"limma"`, `"t_test"`, or `"wilcox_test"`. | `"limma"` |
| `alternative` | Test direction: `"two.sided"`, `"less"`, or `"greater"`. | `"two.sided"` |
| `p_combine_method` | p-value integration method: `"fisher"`, `"z.transform"`, `"logit"`, `"cct"`, `"sumz"`, or `"geometric_mean"`. | `"fisher"` |
| `using_FC` | Logical, used under `alternative = "two.sided"`. If `TRUE`, weight p-values with the actual log2FC; if `FALSE`, weight by `sign(log2FC)`. | `FALSE` |
| `na_ratio` | Maximum allowed missing-value proportion per gene across samples. Genes above this are removed. | `0.5` |
| `using_KNN` | Logical. If `TRUE`, impute missing values with KNN before computing Score2. | `TRUE` |
| `statistics` | Score2 p-value integration: `"arithmetic_mean"`, `"median"`, `"geometric_mean"`, or `"sumz"`. | `"arithmetic_mean"` |
| `purpose` | `"cleaned"` refines de novo gene sets; `"validated"` evaluates external gene sets with fGSEA. | `"cleaned"` |
| `threshold_type` | Thresholding for selecting high-confidence genes: `"quantile"` or `"rank"`. | `"quantile"` |
| `rank_threshold` | Number of top-ranked genes retained when `threshold_type = "rank"`. | `100` |
| `quantile_threshold` | Score2 quantile cutoff (0 to 1). Under `"two.sided"`, genes above and below this quantile are kept. | `0.99` |
| `activation_geneset` | Data frame of activation gene sets to clean (required for `purpose = "cleaned"`). | `NA` |
| `inhibition_geneset` | Data frame of inhibition gene sets to clean (required for `purpose = "cleaned"`). | `NA` |
| `geneSets_gmt` | External gene sets for validation; a data frame with `term` and `gene` columns (for example from `clusterProfiler::read.gmt`). Required for `purpose = "validated"`. | `NA` |
| `min.sz` | Minimum gene set size for fGSEA. | `1` |
| `max.sz` | Maximum gene set size for fGSEA. | `10000` |
| `export_file` | Logical. If `TRUE`, write the cleaned tables to `file_paths`. | `FALSE` |

For `purpose = "cleaned"`, returns
`list(update_activation_geneset, update_inhibition_geneset, cleaning_system)`.
For `purpose = "validated"`, returns `list(fGSEA_result, validated_system)`.

### `Merge_Wnt_genesets()`

Drops too-small gene sets and merges highly similar ones.

```r
Merge_Wnt_genesets(
  file_paths,
  activation_geneset,
  inhibition_geneset,
  delete_GN = TRUE,
  delete_time = "before",
  min_GN = 5,
  integration_method = "Jaccard",
  de_redundant_basis = "igraph_component",
  similarity = 0.3,
  export_file = FALSE
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `file_paths` | Working directory where result files are written when `export_file = TRUE`. | required |
| `activation_geneset` | Data frame of activation gene sets to merge (columns are gene sets). | required |
| `inhibition_geneset` | Data frame of inhibition gene sets to merge. | required |
| `delete_GN` | Logical. If `TRUE`, remove gene sets with too few genes. | `TRUE` |
| `delete_time` | When to filter by size: `"before"` or `"after"` integration. | `"before"` |
| `min_GN` | Minimum number of genes in a gene set; sets below this are removed. | `5` |
| `integration_method` | Similarity coefficient: `"Jaccard"`, `"Sorensen-Dice"`, `"Hub-Promoted"`, or `"Hub-Depressed"`. | `"Jaccard"` |
| `de_redundant_basis` | Redundancy-removal strategy: `"igraph_component"`, or a hierarchical-clustering method such as `"agglomeration_ward.D"`, `"agglomeration_ward.D2"`, `"agglomeration_single"`, `"agglomeration_complete"`, `"agglomeration_average"`, `"agglomeration_mcquitty"`, `"agglomeration_median"`, or `"agglomeration_centroid"`. | `"igraph_component"` |
| `similarity` | Similarity threshold above which gene sets are merged. | `0.3` |
| `export_file` | Logical. If `TRUE`, write the merged tables to `file_paths`. | `FALSE` |

Returns `list(joint_activation_geneset, joint_inhibition_geneset)`; when
`de_redundant_basis = "igraph_component"`, it also returns
`jaccard_network_UP` and `jaccard_network_DN`.

## Notes

1. Data must be accessible in the calling environment and named exactly
   `<Dataset>` and `<Dataset>_G`.
2. Use at least two datasets, as in the example.
3. `file_paths` is only used for writing result files when `export_file =
   TRUE`.

## Author

Dingkang Zhao (赵定康)

Email: <dingkang.25@intl.zju.edu.cn>
