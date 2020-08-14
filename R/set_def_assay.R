#' Set Default Assay
#'
#' This will set and get the default assay to operate on for downstream analysis.
#'
#' @param input the input sce
#' @param assay the character name of the default assay
#' @export
#' @details
#'
#' @examples
#' sce <- set_def_assay(sce, "counts")

set_def_assay <- function(input,
                          def_assay){

  .check_assay(input, def_assay)
  int_metadata(input)$default_assay <- def_assay
  return(input)
}

# Doc Check


