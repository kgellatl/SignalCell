#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input sce object
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' sce <- calc_libsize(input = sce)

calc_libsize <- function(input,
                         assay = "counts",
                         suffix = NULL){

  libsizes <- apply(input@assays@data[[assay]],2,sum)

  if(!is.null(suffix)){
    label <- paste0("libsize", "_", suffix)
  } else {
    label <- "libsize"
  }

  colData(input)[,label] <- libsizes

  return(input)
}
