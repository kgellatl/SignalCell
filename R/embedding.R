#' Get embeddings
#'
#' Returns the 2D matrix associated with a specific embedding
#'
#' @param input the input sce
#' @param the lem within reducedDim()
#' @param the embedding within the LEM
#' @export
#' @details
#'
#' @examples
#' sce <- get_def_assay(sce)

embedding <- function(input,
                      lem,
                      embedding){

  mat <- embeddings(input)

  if(!(lem %in% mat[,"lem"]) || !(embedding %in% mat[,"embedding_key"])){
    stop("Embedding not found, run embeddings() to see available embeddings")
  }

  dims <- metadata(reducedDim(input, lem))$embedding[[embedding]]$embedding
  return(dims)
}

