#' Construct Single Cell Experiment Object
#'
#' This function will take an input sparse expression matrix and make a Single Cell Experiment Object
#'
#' @param input the input data matrix. Provide a sparse matrix
#' @export
#' @details
#' This will take an input matrix and create the Single Cell Experiment Object for further analysis
#' @examples
#' sce <- construct_sce(input = sc_dat)

construct_sce <- function(input) {

  sce <- SingleCellExperiment(list(counts=input))

  sce <- set_def_assay(sce, "counts")

  return(sce)

}

# Doc Check
