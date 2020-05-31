#' Subset Genes
#'
#' This will select genes based on minimum expression, coefficient of variation,
#' or by a preliminary PCA.
#'
#' @param input the input sce
#' @param method can either be "Expression", CV", "PCA", or "Gini"
#' @param assay if NULL will default to def_assay. Can provide a character assay argument.
#' @param threshold UMI threshold for gene detection
#' @param minCells number of cells expressed above threshold for a given gene
#' @param cutoff the percentile of genes to keep
#' @param nComp if method = PCA, the number of components to keep
#' @param log whether or not to log scale the data
#' @param output if "simple" a gene vector will be returned. If "sce" these genes will be annotated to the sce object
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

var_genes <- function(input, method, assay = NULL, threshold = 1, minCells = 10, cutoff = 0.85, nComp = 10, log = F, output = "simple", fudge = F, fudge_val = .01){

  if(is.null(assay)){
    def_assay <- get_def_assay(input)
  }

  input_mat <- input@assays@data@listData[[def_assay]]

  gCount <- apply(input_mat,1,function(x) length(which(x>=threshold))) # a bit wasteful if threshold = 0, but alas.
  gene_subset <- rownames(input_mat[(which(gCount >= minCells)),])

  if(method =="Expression"){

   rowData(input)$var_genes_expression <- NA
   ind <- match(gene_subset, rownames(input))
   rowData(input)$var_genes_expression[ind] <- T

  }

  if(method == "CV"){

    rowData(input)$var_genes_CV <- F

    if(fudge){
      input_mat <- input_mat+fudge_val
    }

    g_exp <- log2(input_mat[gene_subset,]+2)-1
    gsd <- apply(g_exp,1,sd)
    CV <- sqrt((exp(gsd))^2-1)

    ind <- match(names(CV), rownames(input))

    rowData(input)$CV <- NA
    rowData(input)$CV[ind] <- as.vector(CV)

    cv_thresh <- quantile(CV,cutoff)
    gene_subset_cv <- names(which(CV>cv_thresh))

    ind <- match(gene_subset_cv, rownames(input))

    rowData(input)$var_genes_CV[ind] <- T

    gene_subset <- gene_subset_cv

  }

  if(method == "PCA"){
    rowData(input)$var_genes_PCA <- F

    if(fudge){
      input_mat <- input_mat+fudge_val
    }

    if(log){
      input_mat <- log2(input_mat[gene_subset,]+2)-1
    }

    input_scale <- scale(input_mat[gene_subset,])
    pc <- irlba::prcomp_irlba(t(input_scale), nComp, center = F)
    rownames(pc$rotation) <- gene_subset
    d <- mahalanobis(pc$rotation[,1:nComp], center=rep(0, nComp), cov = cov(pc$rotation[,1:nComp]))

    ind <- match(names(d), rownames(input))

    rowData(input)$malhanobis_d <- NA
    rowData(input)$malhanobis_d[ind] <- as.vector(d)

    dThresh <- quantile(d,cutoff)
    gene_subset_malhanobis <- names(which(d>dThresh))

    ind <- match(gene_subset_malhanobis, rownames(input))

    rowData(input)$var_genes_PCA[ind] <- T

    gene_subset <- gene_subset_malhanobis

  }

  if(method == "Gini"){
    rowData(input)$var_genes_gini <- F

    if(fudge){
      input_mat <- input_mat+fudge_val
    }

    if(log){
      input_mat <- log2(input_mat[gene_subset,]+2)-1
    } else {
      input_mat <- input_mat[gene_subset,]
    }

    gini_scores <- edgeR::gini(t(as.matrix(input_mat)))

    ind <- match(names(gini_scores), rownames(input))

    rowData(input)$gini <- NA
    rowData(input)$gini[ind] <- as.vector(gini_scores)

    gini_thresh <- quantile(gini_scores,cutoff)
    gene_subset_gini <- names(which(gini_scores>gini_thresh))

    ind <- match(gene_subset_gini, rownames(input))

    rowData(input)$var_genes_gini[ind] <- T

    gene_subset <- gene_subset_gini


  }

  if(output == "simple"){
    return(gene_subset)
  }
  if(output == "sce"){
    return(input)
  }
}



