#'  embeddings
#'
#' A matrix of available LEM and associated embeddings
#'
#' @param input the input sce
#' @export
#' @details
#'
#' @examples
#' sce <- get_def_assay(sce)

embeddings <- function(input){
    redDims <- names(reducedDims(input))
    embed_mat <- c()
  for (i in 1:length(redDims)) {
    int_dim <- redDims[i]
    embed <- names(metadata(reducedDim(input, int_dim))$embeddings)
    if(!is.null(embed)){
      embed_mat <- c(embed_mat, paste0(int_dim, "-", embed))
    }
  }

    if(!is.null(embed_mat)){
      embed_mat <- matrix(unlist(strsplit(embed_mat, split = "-")), ncol = 2, byrow = T)
      colnames(embed_mat) <- c("lem", "embedding_key")
      rownames(embed_mat) <- paste0("Embedding_", seq(1:nrow(embed_mat)))
      embed_mat <- data.frame(embed_mat)
      rownames(embed_mat) <- paste0("Embedding_", seq(1:nrow(embed_mat)))

    } else {
      embed_mat <- matrix(ncol=2)
      colnames(embed_mat) <- c("lem", "embedding_key")
      embed_mat <- data.frame(embed_mat)
    }

    return(embed_mat)

}

