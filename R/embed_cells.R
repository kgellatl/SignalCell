#' Embed Cells
#'
#' This function will do dimensionality reduction.
#'
#' @param input the input sce
#' @param lem a LEM within reducedDims() must be provided
#' @param embedding_key The name of the embedding within metadata(reducedDim(input)) slot of SCE
#' @param method the method for embedding. tSNE or OTHER?
#' @param comp_select a vector of indices to use for embedding. NULL uses all, seq(1:10) would use first 10 sampleFactors withing the reducedDim(lem)
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended.
#' @param iterations The number of iterations for tSNE to perform.
#' @param write_colData whether or not to write the embedding to colData
#' @param colData_labs a vector of length 2 for the names of the X and Y colnames within colData()
#' @param seed For tSNE, the seed. Can set to NULL if desired.
#' @param print_progress will print progress if TRUE

#' @export
#' @details
#'
#' @examples
#'
embed_cells <- function(input,
                        lem,
                        embedding_key = NULL,
                        method = "tSNE",
                        comp_select = NULL,
                        tSNE_perp = 30,
                        iterations = 500,
                        write_colData = T,
                        colData_labs = NULL,
                        seed = 100,
                        print_progress = T){

  if(!(lem %in% reducedDimNames(input))){
    stop(paste0("LEM not found, LEMs available are ", paste0(reducedDimNames(input), collapse = ", ")))
  }

  if(is.null(embedding_key)){
    embedding_key <- method
  }

  if(is.null(colData_labs)){
    colData_labs <- paste0(method, c("_x", "_y"))
  }

  input_lem <- reducedDim(input, lem)
  embed_met <- metadata(input_lem)
  tsne_input <- sampleFactors(input_lem)

  if(!is.null(comp_select)){
    if(max(comp_select) > ncol(tsne_input)){
      print(paste0("lem embedding only contains "), ncol(tsne_input), " dimensions. max(comp_select) must be less than this value.")
    }
    tsne_input <- tsne_input[,comp_select]
  } else {
    comp_select <- seq(1:ncol(tsne_input))
  }

  if(method == "tSNE"){
    if(print_progress == TRUE){
      print("Starting tSNE")
    }

    set.seed(seed)
    embedding <- Rtsne::Rtsne(tsne_input, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
    set.seed(NULL)
    embedding <- embedding$Y
    row.names(embedding) <- rownames(tsne_input)
    colnames(embedding) <- c("x", "y")
    embedding[,"x"] <-  abs(min(embedding[,"x"]))+embedding[,"x"]
    embedding[,"x"] <-  embedding[,"x"]/max(embedding[,"x"])
    embedding[,"y"] <-  abs(min(embedding[,"y"]))+embedding[,"y"]
    embedding[,"y"] <-  embedding[,"y"]/max(embedding[,"y"])

    colnames(embedding) <- colData_labs
  }

  if(write_colData){
    colData(input)[colData_labs[1]] <- embedding[,colData_labs[1]]
    colData(input)[colData_labs[2]] <- embedding[,colData_labs[2]]
  }

  ####

  parameters <- list(lem, method, embedding_key, comp_select, tSNE_perp, iterations, write_colData, colData_labs, seed)
  names(parameters) <- c("lem", "method", "embedding_key", "comp_select", "tSNE_perp", "iterations", "write_colDat", "colDat_labs", "seed")

  embed_metadat <- list(embedding, parameters)
  names(embed_metadat) <- c("embedding", "parameters")
  embed_metadat <- list(embed_metadat)
  names(embed_metadat) <- embedding_key

  if(embedding_key %in% embeddings(input)$embedding_key){
    metadata(reducedDim(input, lem))[["embeddings"]][embedding_key] <- embed_metadat
  } else {
    metadata(reducedDim(input, lem))[["embeddings"]] <- append(metadata(reducedDim(input, lem))[["embeddings"]], embed_metadat)
  }

  return(input)

}
