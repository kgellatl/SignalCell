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
#' @param write_colDat whether or not to write the embedding to colData
#' @param colDat_labs a vector of length 2 for the names of the X and Y colnames within colData()
#' @param seed For tSNE, the seed. Can set to NULL if desired.
#' @param print_progress will print progress if TRUE

#' @export
#' @details
#'
#' @examples
#'
embed_cells <- function(input,
                        lem,
                        embedding_key,
                        method = "tSNE",
                        comp_select = NULL,
                        tSNE_perp = 30,
                        iterations = 500,
                        write_colDat = T,
                        colDat_labs = c("x","y"),
                        seed = 100,
                        print_progress = T){

  if(!(lem %in% reducedDimNames(input))){
    stop(paste0("LEM not found, LEMs available are ", paste0(reducedDimNames(input), collapse = ", ")))
  }
  input_lem <- reducedDim(input, lem)
  embed_met <- metadata(input_lem)
  tsne_input <- sampleFactors(input_lem)

  if(!is.null(comp_select)){
    if(max(comp_select) > ncol(tsne_input)){
      print(paste0("lem embedding only contains "), ncol(tsne_input), " dimensions. max(comp_select) must be less than this value")
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
  }

  if(write_colDat){
    colData(input)[colDat_labs[1]] <- embedding[,"x"]
    colData(input)[colDat_labs[2]] <- embedding[,"y"]
  }


  ####

  parameters <- list(lem, method, embedding_key, comp_select, tSNE_perp, iterations, write_colDat, colDat_labs, seed)
  names(parameters) <- c("lem", "method", "embedding_key", "comp_select", "tSNE_perp", "iterations", "write_colDat", "colDat_labs", "seed")

  embed_metadat <- list(embedding, parameters)
  names(embed_metadat) <- c("embedding", "parameters")
  embed_metadat <- list(embed_metadat)
  names(embed_metadat) <- embedding_key

  ### Combining with existing LEM metadata
  str(embed_metadat) # the embedding
  str(embed_met) # the lem


  embed_met <- append(embed_met, embed_metadat)
  str(embed_met)

  return(input)

}
