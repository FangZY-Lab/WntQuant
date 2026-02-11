# WntQuant

## Overview

**WntQuant: Robust Quantification of Wnt/β-catenin Pathway Activity through Multi-Dataset Integration: A Context-Insensitive and Directional Approach for Precision Oncology**

WntQuant is an R package that enables systematic derivation and refinement of Wnt pathway-associated gene signatures from heterogeneous transcriptomic data. The framework implements a dual-phase strategy combining cross-study differential meta-analysis with evidence-based purification to establish high-confidence directional Wnt pathway signatures.

**Get_Wnt_denovo_genesets: De Novo Wnt Signature Derivation Module**  
A comprehensive meta-analytic pipeline that identifies Wnt pathway-associated genes through:
- **Multi-dataset differential analysis**: Implements limma, t-test, and Wilcoxon tests across independent bulk RNA-seq studies stratified by Wnt/β-catenin activity
- **Directional evidence synthesis**: Designates up-regulated and down-regulated genes as activation and inhibition signatures respectively
- **Flexible p-value integration**: Supports Fisher, z-transform, logit, Cauchy combination test (CCT), sumz, and geometric mean methods for cross-study evidence aggregation
- **Adaptive gene selection**: Offers rank-based (top N genes) or threshold-based (p-value + log2FC) strategies for signature derivation
- **Scalable architecture**: Processes multiple independent datasets with automatic grouping and result consolidation

**Wnt_purification_system: Precision Signature Refinement Module**  
A robust validation and cleaning system that enhances Wnt signature accuracy through:
- **Score2 metric computation**: Integrates directional p-values and fold changes across studies into a quantitative confidence score
- **Intelligent missing value handling**: Optional KNN imputation with configurable missing rate thresholds
- **Dual-purpose functionality**: Supports both de novo signature cleaning and external Wnt gene set validation
- **Evidence-based prioritization**: Employs quantile or rank-based thresholds to retain high-confidence Wnt-associated genes
- **fGSEA integration**: Enables systematic evaluation of public Wnt gene sets against user-defined expression compendia

**Merge_Wnt_genesets: Biological Redundancy Resolution Module**  
A network-based integration tool that consolidates biologically similar Wnt signatures through:
- **Multiple similarity metrics**: Implements Jaccard, Sørensen-Dice, Hub-Promoted, and Hub-Depressed coefficients for gene set comparison
- **Flexible clustering strategies**: Supports graph component detection and eight agglomerative hierarchical clustering methods (ward.D, ward.D2, single, complete, average, mcquitty, median, centroid)
- **Adaptive size filtering**: Optional removal of gene sets below minimum gene count thresholds (pre- or post-integration)
- **Union-based consolidation**: Merges highly similar gene sets through intersection-weighted union operations
- **Network visualization ready**: Outputs similarity matrices for downstream Jaccard network construction

**Framework Features:**
- De novo Wnt signature discovery from multi-study expression compendia
- Quantitative confidence scoring (Score2) for gene-level Wnt association
- Context-specific signature refinement across diverse biological settings
- Comprehensive support for one-sided and two-sided testing frameworks
- Seamless integration with downstream GSEA and fGSEA workflows
- Reproducible outputs in standardized tabular formats

WntQuant provides a systematic and statistically rigorous approach for establishing directional Wnt pathway signatures, offering improved specificity and biological interpretability compared to single-dataset or undirected enrichment methods.

## Installation

### Prerequisites

Make sure you have R (version 4.0.0 or higher) installed. The package requires several Bioconductor and CRAN dependencies.

### Installation Steps

```r
# Install from GitHub
devtools::install_github("FangZY-Lab/WntQuant", force = TRUE)
```

## Quick Start

### De Novo Wnt Signature Discovery

```r
# Identify Wnt activation and inhibition signatures from multiple datasets
denovo_sigs <- Get_Wnt_denovo_genesets(
  file_paths = "./data",
  expression_accession_vector = c("GSEA_breast", "GSEB_colorectal"),
  group_HL = wnt_group_assignments,
  gene_difference_method = "limma",
  alternative = "two.sided",
  p_combine_method = "cct",
  threshold_or_rank = "rank",
  top_genes = 500
)

# Access results
wnt_activation <- denovo_sigs$activation
wnt_inhibition <- denovo_sigs$inhibition
```

### Signature Purification and Validation

```r
# Refine derived signatures with Score2 metric
refined_sigs <- Wnt_purification_system(
  file_paths = "./data",
  expression_accession_vector = c("GSEA_breast", "GSEB_colorectal"),
  group_HL = wnt_group_assignments,
  alternative = "two.sided",
  purpose = "cleaned",
  threshold_type = "quantile",
  quantile_threshold = 0.99,
  activation_geneset = denovo_sigs$activation,
  inhibition_geneset = denovo_sigs$inhibition
)

# Validate external Wnt gene sets
validation_result <- Wnt_purification_system(
  file_paths = "./data",
  expression_accession_vector = c("GSEA_breast", "GSEB_colorectal"),
  group_HL = wnt_group_assignments,
  alternative = "two.sided",
  purpose = "validated",
  geneSets_gmt = "wnt_pathways.gmt"
)
```

### Merge Redundant Signatures

```r
# Consolidate biologically similar gene sets
merged_sigs <- Merge_Wnt_genesets(
  file_paths = "./results",
  activation_geneset = refined_sigs$update_activation_geneset,
  inhibition_geneset = refined_sigs$update_inhibition_geneset,
  delete_GN = TRUE,
  min_GN = 5,
  integration_method = "Jaccard",
  de_redundant_basis = "igraph_component",
  similarity = 0.3
)

# Access merged signatures
final_activation <- merged_sigs$joint_activation_geneset
final_inhibition <- merged_sigs$joint_inhibition_geneset
jaccard_network <- merged_sigs$jaccard_network_UP
```

## Input Data Format

### Expression Data
- Format: Data frames or matrices with genes in rows, samples in columns
- Row names: Gene identifiers (e.g., ENTREZID, SYMBOL)
- Column names: Sample identifiers matching group annotation files

### Group Annotation (`group_HL`)
A data frame with the following columns:
- `Accession`: Dataset identifier matching expression_accession_vector
- `Group1`: Name of first comparison group
- `Group2`: Name of second comparison group  
- `Group1_Status`: Wnt activity status of Group1 ("H" or "L")
- `Group2_Status`: Wnt activity status of Group2 ("H" or "L")

### Sample Group Files (`[Accession]_G`)
For each dataset, a corresponding group file named `[Accession]_G` containing:
- `Tag`: Sample identifiers matching expression data column names
- `group`: Group assignment matching Group1/Group2 in group_HL

## License

This package is licensed under the GPL-3.0 License.

## Contact

For questions, bug reports, or feature requests, please contact the corresponding author(dingkang.22@intl.zju.edu.cn) or open an issue on GitHub.
