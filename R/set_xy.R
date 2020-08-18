#' Set xy
#'
#' This will set default X and Y coordinates based on the LEM results within reducedDim(sce)
#'
#' @param input the input sce
#' @param lem the LEM within reducedDims
#' @param x The x values for 2D representation. Can be character x or a numeric column of sampleFactors(lem)
#' @param y The y values for 2D representation. Can be character y or a numeric column of sampleFactors(lem)
#' @export
#' @details
#'
#' @examples
#' sce <- set_xy(input = sce, lem = "iPCA", x = "x", y = "y")

set_xy <- function(input,
                   lem,
                   embedding,
                   x = "x",
                   y = "y"){

  if(!(lem %in% names(reducedDims(input)))){
    stop(paste0("LEM not found, reducedDims available are, ", names(reducedDims(input))))
  }

  if(!(embedding %in% embeddings(input)$embedding_key)){
    stop(paste0("Embedding not found, used embeddings() to view available options."))
  }

  embed <- embedding(input, lem, embedding)
  colnames(embed) <- c(x,y)

  colData(input)[,x] <- embed[,x]
  colData(input)[,y] <- embed[,y]
  return(input)
}

# Doc Check


