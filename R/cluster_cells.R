#' Cluster Single Cell
#'
#' This will perform clustering on your single cell data.
#'
#' @param input the input ex_sc
#' @param lem the LEM within ReducedDims slot
#' @param dims either "2d" or "Comp"
#' @param method can either be "spectral" or "density" which is on 2d
#' @param num_clust the number of clusters. Required for spectral but optional for density.
#' @param name name of the colData cluster column
#' @param s the number of standard deviations from the curve to select cluster centers
#' @export
#' @details
#' This will perform clustering on either the high dimensional PCA / ICA components if dimension = Comp,
#' or the 2d tsne result if method = density. Typically spectral clustering works much better on higher dimensional data,
#' which density based clustering works better on 2d data.
#' @examples
#' ex_sc_example <- cluster_sc(input = ex_sc_example, dimension = "Comp", method = "spectral", num_clust = 6)

cluster_cells <- function(input,
                          lem,
                          dims,
                          method,
                          num_clust = NULL,
                          name = "Cluster",
                          s=2) {

  if(!(lem %in% reducedDimNames(sce))){
    stop(paste0("LEM not found, LEMs available are ", paste0(reducedDimNames(input), collapse = ", ")))
  }

  if(!(dims %in% c("2d", "Comp"))){
    stop("dims must be either '2d' or 'Comp'")
  }

  if(!(method %in% c("spectral", "density"))){
    stop("method must be either 'spectral' or 'density'")
  }

  dim_dat <- reducedDim(input, lem)

  if(dims == "2d"){
    x <- dim_dat@metadata$x
    y <- dim_dat@metadata$y
    tocluster <- matrix(c(x,y), ncol = 2)
    rownames(tocluster) <- colnames(input)
  } else {
    tocluster <- dim_dat@sampleFactors
  }

  if(method == "spectral"){
    if(is.null(num_clust)){
      stop("Please provide num_clust argument")
    }
    spec <- kknn::specClust(tocluster, centers = num_clust, method = 'random-walk')
    sc_clusters <- spec$cluster
  }

  if(method == "density"){
    dist_tSNE = dist(tocluster)
    clust_tSNE = densityClust::densityClust(dist_tSNE, gaussian = T) # calculate density
    comb = as.data.frame(clust_tSNE$rho*clust_tSNE$delta) # combine rho and delta values
    colnames(comb) = "gamma"
    comb = comb[order(comb$gamma, decreasing = T), ,drop=F]
    comb$index = seq(nrow(comb))
    if (is.null(num_clust)) {
      # chose the max k from gamma distribution
      fit = mgcv::gam(formula = gamma ~ s(index, bs="cs"), data = log10(comb[floor(0.01*nrow(comb)):nrow(comb),]+1))
      test = log10(comb[,"index", drop=F])
      p = predict(fit, test, type = "link", se.fit = T)
      comb$pred = (10^predict(fit, test))-1
      comb$residual = comb$gamma-comb$pred
      comb$predsd = comb$pred+(s*sd(comb$residual))
      print(ggplot(comb, aes(index, gamma)) +
              geom_point(size = 0.5) +
              theme_bw() +
              scale_x_log10() +
              scale_y_log10() +
              geom_line(data=comb, aes(index, pred), colour="red") +
              geom_line(data=comb, aes(index, predsd), colour="blue"))
      cellcut = rownames(comb[comb$gamma>comb$predsd,])
    } else {
      cellcut = rownames(comb)[1:num_clust]
    }
    cellidx = which(names(clust_tSNE$rho)%in%cellcut)
    sc_clusters = densityClust::findClusters(clust_tSNE, peaks = cellidx)
    sc_clusters <- sc_clusters$clusters
  }

  colData(input)[name] <- paste0("Cluster_", sc_clusters)
  return(input)

}
