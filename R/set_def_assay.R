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

  if(!(def_assay %in% assayNames(input))){
    stop(paste0("Assay not found, cannot set to default, assays available are, ", names(input@assays)))
  }
  input@assays@data@metadata$default_assay <- def_assay
  return(input)
}

# Doc Check


