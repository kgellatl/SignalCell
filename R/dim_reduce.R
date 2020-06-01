#' Dimension Reduction
#'
#' This function will do dimensionality reduction.
#'
#' @param input the input sce
#' @param genelist the subset of genes to perform dimensionality reduction on
#' @param assay The assay to operate on. Will default to get_def_assay(input)
#' @param pre_reduce the algorithm choice for reduction before tSNE (either "ICA", "PCA", "iPCA").
#' @param nComp the number of components to reduce too before tSNE, 5-20 recommended.
#' @param tSNE_perp number of cells expressed above threshold for a given gene, 10-100 recommended.
#' @param iterations The number of iterations for tSNE to perform.
#' @param print_progress will print progress if TRUE
#' @param nVar cutoff for percent of variance explained from PCs
#' @param log Whether or not to log  the input assay
#' @param scale Whether or not to scale the input assay
#' @param save_lem Whether or not to write result to ReducedDim slot of SCE
#' @param reducedDim_key The name of the LEM within ReducedDim slot of SCE
#' @param seed For tSNE, the seed. Can set to NULL if desired.
#' @importFrom fastICA fastICA
#' @importFrom  Rtsne Rtsne
#' @importFrom irlba prcomp_irlba
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = gene_subset, pre_reduce = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
dim_reduce <- function(input,
                        genelist,
                        assay = NULL,
                        pre_reduce = "iPCA",
                        nComp = 15,
                        tSNE_perp = 30,
                        iterations = 1000,
                        print_progress=TRUE,
                        nVar=.85,
                        log = F,
                        scale = T,
                        save_lem = T,
                        reducedDim_key = NULL,
                        seed = 100){

  if(is.null(assay)){
    def_assay <- get_def_assay(input)
    input_mat <- input@assays@data@listData[[def_assay]]
    assay_name <- def_assay
  } else {
    if(!(assay%in%names(input@assays))){
      stop(paste0("Assay not found, assays available are, ", names(input@assays)))
    }
    input_mat <- input@assays@data@listData[[assay]]
    assay_name <- assay
  }

  args_list <- list(assay_name, genelist, pre_reduce, nComp, tSNE_perp, iterations, nVar, log, scale)
  names(args_list) <- c("assay_name", "genelist", "pre_reduce", "nComp", "tSNE_perp", "iterations", "nVar", "log", "scale")
  metadata_lem <- args_list

  input_mat <- input_mat[gene_subset,]
  if(log){
    input_mat <- log2(input_mat[,]+2)-1
  }
  if(scale){
    input_mat <- scale(input_mat)
  }

  if(pre_reduce == "ICA"){
    if(print_progress == TRUE){
      print("Starting ICA")
    }
    ica <- fastICA::fastICA(t(input_mat), n.comp=(nComp), alg.typ = 'parallel', fun='logcosh', alpha = 1.0, method = 'C', verbose = print_progress)
    colnames(ica$A) <- gene_subset
    rownames(ica$S) <- colnames(input)
    colnames(ica$S) <- paste0("IC_Comp", seq(1:ncol(ica$S)))

    tsne_input = ica$S

    sampleFactors_lem <- ica$S
    featureLoadings_lem <- t(ica$A)
    factorData_lem <- DataFrame(ica$W)
  }

  if(pre_reduce == "PCA"){
    if(print_progress == TRUE){
      print("Starting PCA")
    }
    PCA <- irlba::prcomp_irlba(t(input_mat), nComp, center = F)
    rownames(PCA$x) <- colnames(input)
    colnames(PCA$x) <- paste0("PC_Comp", seq(1:ncol(PCA$x)))

    tsne_input = PCA$x

    sampleFactors_lem <- PCA$x
    featureLoadings_lem <- PCA$rotation
    factorData_lem <- DataFrame(PCA$sdev)
  }

  if(pre_reduce == "iPCA"){
    if(print_progress == TRUE){
      print("Starting iPCA")
    }
    iPCA <- irlba::prcomp_irlba(input_mat, nComp, center = F)
    rownames(iPCA$rotation) <- colnames(input)
    colnames(iPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(iPCA$rotation)))

    tsne_input = iPCA$rotation

    sampleFactors_lem <- iPCA$rotation
    featureLoadings_lem <- iPCA$x
    factorData_lem <- DataFrame(iPCA$sdev)
  }

  if(pre_reduce == "vPCA"){
    if(print_progress == TRUE){
      print("Starting vPCA")
    }
    vPCA <- irlba::prcomp_irlba(input_mat, n = nComp)
    # sum components until variance is >= x%
    var = vPCA$sdev^2/sum(vPCA$sdev^2)
    totalvar = var[1]
    maxPC = 1
    while (totalvar < nVar) {
      maxPC=maxPC+1
      totalvar = sum(var[1:maxPC])
    }

    if(maxPC < 2){
      stop("Percent variance threshold has left than 2 PCs. Please increase this value.")
    }
    vPCA$rotation = vPCA$rotation[,1:maxPC]
    rownames(vPCA$rotation) <- colnames(input)
    colnames(vPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(vPCA$rotation)))

    tsne_input = vPCA$rotation

    sampleFactors_lem <- vPCA$rotation
    featureLoadings_lem <- vPCA$x[,1:maxPC]
    factorData_lem <- DataFrame(vPCA$sdev[1:maxPC])
  }

  if(print_progress == TRUE){
    print("Starting tSNE")
  }

  set.seed(seed)
  tSNE_result <- Rtsne::Rtsne(tsne_input, dims = 2, perplexity = tSNE_perp, theta = 0.5, check_duplicates = F, pca = F, max_iter = iterations, verbose = print_progress)
  set.seed(NULL)
  tSNE_result <- tSNE_result$Y

  row.names(tSNE_result) <- rownames(tsne_input)
  colnames(tSNE_result) <- c("x", "y")
  tSNE_result[,"x"] <-  abs(min(tSNE_result[,"x"]))+tSNE_result[,"x"]
  tSNE_result[,"x"] <-  tSNE_result[,"x"]/max(tSNE_result[,"x"])
  tSNE_result[,"y"] <-  abs(min(tSNE_result[,"y"]))+tSNE_result[,"y"]
  tSNE_result[,"y"] <-  tSNE_result[,"y"]/max(tSNE_result[,"y"])


  colData(input)$x <- tSNE_result[,"x"]
  colData(input)$y <- tSNE_result[,"y"]


  if(save_lem){

    if(is.null(reducedDim_key)){
      reducedDim_key <- pre_reduce
    }

    rownames(featureLoadings_lem) <- genelist

    lem <- LinearEmbeddingMatrix(sampleFactors_lem,
                                 featureLoadings_lem,
                                 factorData_lem,
                                 metadata_lem)

    lem@metadata$x <- colData(input)$x
    lem@metadata$y <- colData(input)$y

    reducedDim(input, type = reducedDim_key) <- lem

  }

  return(input)
}

