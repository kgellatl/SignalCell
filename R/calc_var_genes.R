#' Subset Genes
#'
#' This will select genes based on minimum expression, coefficient of variation,
#' or by a preliminary PCA.
#'
#' @param input the input sce
#' @param method can either be "CV", "Malhanobis", or "Gini"
#' @param assay if NULL will default to def_assay. Can provide a character assay argument.
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @param nComp if method = PCA, the number of components to keep
#' @param log whether or not to log scale the data
#' @param fudge whether or not to add a fudge factor to the entire matrix
#' @param fudge_val the value to add to matrix before transformation
#' @export
#' @details
#' Genes will be first filtered by minimum expression selecting by subsetting to genes that are
#' expressed above the threshold in more than minCells.
#' If the method is CV, it will first subset the genes based on the expression cutoffs,
#' then find the coefficient of variation across all genes.
#' Next it will select the percentile of genes (cutoff) based on their coefficient of variation.
#' The last method will perform PCA on the cells, and then look at the loadings of each gene.
#' By finding genes that are off center (via malhanoobis distance) we can filter to include only genes that contribute significant variance to the data.
#' @examples
#' gene_subset <- subset_genes(input = sce, method = "PCA", assay = "counts")

calc_var_genes <- function(input, method, assay = NULL, threshold = 1, minCells = 10, nComp = 10, log = F, fudge = F, fudge_val = .01){

  if(is.null(assay)){
    def_assay <- get_def_assay(input)
  }

  if(!(method%in%c("CV", "Malhanobis", "Gini"))){
    stop("Method must be one of CV, Malhanobis, or Gini")
  }
  input_mat <- input@assays@data@listData[[def_assay]]

  gCount <- apply(input_mat,1,function(x) length(which(x>=threshold)))
  gene_subset <- rownames(input_mat[(which(gCount >= minCells)),])


  if(method == "CV"){
    if(fudge){
      input_mat <- input_mat[gene_subset,]+fudge_val
    }

    g_exp <- log2(input_mat[gene_subset,]+2)-1
    gsd <- apply(g_exp,1,sd)
    CV <- sqrt((exp(gsd))^2-1)

    ind <- match(names(CV), rownames(input))
    rowData(input)$CV <- NA
    rowData(input)$CV[ind] <- as.vector(CV)

  }

  if(method == "Malhanobis"){
    if(fudge){
      input_mat <- input_mat[gene_subset,]+fudge_val
    }
    if(log){
      input_mat <- log2(input_mat[gene_subset,]+2)-1
    }

    input_scale <- scale(input_mat[gene_subset,])
    pc <- irlba::prcomp_irlba(t(input_scale), nComp, center = F)
    rownames(pc$rotation) <- gene_subset
    d <- mahalanobis(pc$rotation[,1:nComp], center=rep(0, nComp), cov = cov(pc$rotation[,1:nComp]))

    ind <- match(names(d), rownames(input))
    rowData(input)$Malhanobis <- NA
    rowData(input)$Malhanobis[ind] <- as.vector(d)

  }

  if(method == "Gini"){
    if(fudge){
      input_mat <- input_mat[gene_subset,]+fudge_val
    }
    if(log){
      input_mat <- log2(input_mat[gene_subset,]+2)-1
    }

    gini_scores <- edgeR::gini(t(as.matrix(input_mat[gene_subset,])))

    ind <- match(names(gini_scores), rownames(input))
    rowData(input)$Gini <- NA
    rowData(input)$Gini[ind] <- as.vector(gini_scores)
  }

  return(input)
}



