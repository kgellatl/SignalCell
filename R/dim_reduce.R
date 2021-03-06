#' Dimension Reduction
#'
#' This function will do dimensionality reduction.
#'
#' @param input the input sce
#' @param genelist the subset of genes to perform dimensionality reduction on
#' @param method the algorithm choice for reduction before tSNE (either "ICA", "PCA", "iPCA", or FALSE if you want to reuse).
#' @param assay The assay to operate on. Will default to get_def_assay(input)
#' @param nComp the number of components to reduce too before tSNE, 5-20 recommended.
#' @param nVar cutoff for percent of variance explained from PCs
#' @param log Whether or not to log  the input assay
#' @param scale Whether or not to scale the input assay
#' @param reducedDim_key The name of the LEM within ReducedDim slot of SCE
#' @param seed For tSNE, the seed. Can set to NULL if desired.
#' @importFrom fastICA fastICA
#' @importFrom  Rtsne Rtsne
#' @importFrom irlba prcomp_irlba
#' @export
#' @details
#' If the method is ICA, independent component analysis will be performed, and then tSNE will do the final dimension reduction. If PCA is selected, PCA will be performed before on the expression matrix transpose before tSNE. This PCA will use the cells positions on the principal components. If iPCA is selected, PCA will be be performed but without transposing the data. This will create "meta cells" instead of meta genes created in the typical PCA. Then tSNE will be performed on each cells contribution (loading) to the meta cell. We find that iPCA is much more robust and leads to cleaner clusters than traditional PCA.
#' @examples
#' ex_sc_example <- dim_reduce(input = ex_sc_example, genelist = genelist, method = "iPCA", nComp = 15, tSNE_perp = 30, iterations = 500, print_progress=TRUE)
#'
dim_reduce <- function(input,
                       genelist,
                       method,
                       assay = NULL,
                       nComp = 15,
                       nVar=.85,
                       log = F,
                       scale = T,
                       reducedDim_key = NULL,
                       seed = 100){

  if(is.null(assay)){
    def_assay <- get_def_assay(input)
    input_mat <- assay(input, def_assay)[genelist,]
    assay_name <- def_assay
  } else {
    .check_assay(input, assay)
    input_mat <- assay(input, assay)[genelist,]
    assay_name <- assay
  }

  if(!(method%in%c("ICA", "PCA", "iPCA", "vPCA", FALSE))){
    stop("Pre-reduce muse be one of ICA, PCA, iPCA, vPCA, or FALSE")
  }

  args_list <- list(assay_name, genelist, method, nComp, nVar, log, scale)
  names(args_list) <- c("assay_name", "genelist", "method", "nComp", "nVar", "log", "scale")
  args_list <- list(args_list)
  names(args_list) <- "parameters"
  metadata_lem <- list(args_list, list())
  names(metadata_lem) <- c("lem", "embeddings")

  if(log){
    input_mat <- log2(input_mat+2)-1
  }
  if(scale){
    input_mat <- scale(input_mat)
  }

  if(method == "ICA"){
    ica <- fastICA::fastICA(t(input_mat), n.comp=(nComp), alg.typ = 'parallel', fun='logcosh', alpha = 1.0, method = 'C', verbose = print_progress)
    colnames(ica$A) <- genelist
    rownames(ica$S) <- colnames(input)
    colnames(ica$S) <- paste0("IC_Comp", seq(1:ncol(ica$S)))

    sampleFactors_lem <- ica$S
    featureLoadings_lem <- t(ica$A)
    factorData_lem <- DataFrame(ica$W)
  }

  if(method == "PCA"){
    PCA <- irlba::prcomp_irlba(t(input_mat), nComp)
    rownames(PCA$x) <- colnames(input)
    colnames(PCA$x) <- paste0("PC_Comp", seq(1:ncol(PCA$x)))

    sampleFactors_lem <- PCA$x
    featureLoadings_lem <- PCA$rotation
    factorData_lem <- DataFrame(PCA$sdev)
  }

  if(method == "iPCA"){
    iPCA <- irlba::prcomp_irlba(input_mat, nComp)
    rownames(iPCA$rotation) <- colnames(input)
    colnames(iPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(iPCA$rotation)))

    sampleFactors_lem <- iPCA$rotation
    featureLoadings_lem <- iPCA$x
    factorData_lem <- DataFrame(iPCA$sdev)
  }

  if(method == "vPCA"){
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
      stop("Percent variance threshold has left less than 2 PCs. Please increase this value.")
    }

    vPCA$rotation = vPCA$rotation[,1:maxPC]
    rownames(vPCA$rotation) <- colnames(input)
    colnames(vPCA$rotation) <- paste0("iPC_Comp", seq(1:ncol(vPCA$rotation)))

    sampleFactors_lem <- vPCA$rotation
    featureLoadings_lem <- vPCA$x[,1:maxPC]
    factorData_lem <- DataFrame(vPCA$sdev[1:maxPC])
  }

  if(is.null(reducedDim_key)){
    reducedDim_key <- method
  }

  rownames(featureLoadings_lem) <- genelist
  lem <- LinearEmbeddingMatrix(sampleFactors_lem,
                               featureLoadings_lem,
                               factorData_lem,
                               metadata_lem)
  reducedDim(input, type = reducedDim_key) <- lem

  return(input)
}




