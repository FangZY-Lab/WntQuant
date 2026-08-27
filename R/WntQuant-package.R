#' @importFrom doBy summaryBy
#' @importFrom dplyr mutate
#' @importFrom fgsea fgsea
#' @importFrom igraph components graph_from_data_frame
#' @importFrom impute impute.knn
#' @importFrom limma contrasts.fit eBayes lmFit makeContrasts topTable
#' @importFrom magrittr "%>%"
#' @importFrom metap sumz
#' @importFrom stats as.dist cutree hclust median model.matrix na.omit pcauchy pt quantile t.test wilcox.test
#' @importFrom survcomp combine.test
#' @importFrom utils write.table
NULL

utils::globalVariables(c("mean1", "mean2", "H", "L", "logFC", "P.Value", "."))
