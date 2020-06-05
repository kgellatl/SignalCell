#' Get Default Assay
#'
#' This will get the default assay used for downstream analysis.
#'
#' @param input the input sce
#' @export
#' @details
#'
#' @examples
#' sce <- get_def_assay(sce)

get_def_assay <- function(input){
  def_assay <- int_metadata(sce)$default_assay
  return(def_assay)
}

# Doc check
