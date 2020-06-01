#' Get Variable Genes
#'
#' This will get the variable genes above a threshold, based on quantile.
#'
#' @param input the input sce
#' @param method One of "CV", "Malhanobis", or "Gini
#' @param cutoff The threshold below which a gene will not be returned. Expressed as decimal.
#' @export
#' @details
#'
#' @examples
#' sce <- get_def_assay(sce)

get_var_genes <- function(input,
                          method,
                          cutoff){

  if(!(method%in%c("CV", "Malhanobis", "Gini"))){
    stop("Method must be one of CV, Malhanobis, or Gini")
  }

  if(!(method%in%colnames(rowData(sce)))){
    stop("Please use calc_var_genes to calculate genewise statistics before this function.")
  }

  vals <- rowData(sce)[method]
  vals <- vals[,1]
  cut_val <- quantile(vals, cutoff, na.rm = T)
  ind <- which(vals >= cut_val)
  gene_subset <- rownames(sce)[ind]

  return(gene_subset)
}
