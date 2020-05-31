#' Get Default Assay
#'
#' This will get the default assay to operate on for downstream analysis.
#'
#' @param input the input sce
#' @export
#' @details
#'
#' @examples
#' sce <- get_def_assay(sce)

get_def_assay <- function(input){
  def_assay <- input@assays@data@metadata$default_assay
  return(def_assay)
}
