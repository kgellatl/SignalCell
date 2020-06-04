#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input sce object
#' @param assay which assay to calculate library size on
#' @param suffix the
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' sce <- calc_libsize(input = sce)

calc_libsize <- function(input,
                         assay = "counts",
                         label = NULL){

  if(!(assay %in% assayNames(input))){
    stop(paste0("Assay not found, cannot set to default, assays available are, ", names(input@assays)))
  }

  libsizes <- apply(assay(input, assay),2,sum)

  if(is.null(label)){
    label <- "libsize"
  }

  colData(input)[,label] <- libsizes

  return(input)
}
