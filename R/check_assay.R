#' check_assay
#'
#' This function will check the assay provided
#'
#' @param input the input sce object
#' @param assay which assay to calculate library size on
#' @export
#' @details
#' This will ensure the provided assay is found within assays()
#' @examples
#' .check_assay(input = sce, assay = "counts)

.check_assay <- function(input,
                         assay){

  if(!(assay %in% assayNames(input))){
    stop(paste0("Assay not found. Assays available are : ", paste0(assayNames(input), collapse = " ")))
  }
}
