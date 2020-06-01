#' Set xy
#'
#' This will set default X and Y coordinates based on the LEM results within reducedDims(sce)
#'
#' @param input the input sce
#' @param lem the LEM within reducedDims
#' @param x The x values for 2D representation. Can be character x or a numeric column of lem@sampleFactors
#' @param y The y values for 2D representation. Can be character y or a numeric column of lem@sampleFactors
#' @export
#' @details
#'
#' @examples
#' sce <- set_xy(input = sce, lem = "iPCA", x = "x", y = "y")

set_xy <- function(input,
                   lem,
                   x = "x",
                   y = "y"){

  if(!(lem %in% names(reducedDims(input)))){
    stop(paste0("LEM not found, reducedDims available are, ", names(reducedDims(input))))
  }

  int_lem <- reducedDim(input, lem)

  if(x == "x"){
    colData(input)$x <- int_lem@metadata$x
  } else {
    colData(input)$x <- int_lem@sampleFactors[,x]
  }

  if(y == "y"){
    colData(input)$y <- int_lem@metadata$y
  } else {
    colData(input)$y <- int_lem@sampleFactors[,y]

  }
    return(input)
}



