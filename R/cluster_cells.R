#' Cluster Single Cell
#'
#' This will perform clustering on your single cell data.
#'
#' @param input the input ex_sc
#' @param lem the LEM within ReducedDims()
#' @param dims either "2d" or "Comp"
#' @param method can either be "spectral" or "density" which is on 2d
#' @param embedding if method is 2D, an embedding must be provided
#' @param k the number of clusters. Required for spectral but optional for density.
#' @param name name of the colData cluster column
#' @param s the number of standard deviations from the curve to select cluster centers
#' @export
#' @details
#' This will perform clustering on either the high dimensional PCA / ICA components if dimension = Comp,
#' or the 2d tsne result if method = density. Typically spectral clustering works much better on higher dimensional data,
#' which density based clustering works better on 2d data.
#' @examples
#' ex_sc_example <- cluster_sc(input = ex_sc_example, dimension = "Comp", method = "spectral", k = 6)

cluster_cells <- function(input,
                          lem,
                          dims,
                          method,
                          embedding = NULL,
                          xy = NULL,
                          k = NULL,
                          name = "Cluster",
                          s=2,
                          cluster_stats = F) {

  if(!(dims %in% c("2d", "Comp"))){
    stop("dims must be either '2d' or 'Comp'")
  }

  if(!(method %in% c("spectral", "density"))){
    stop("method must be either 'spectral' or 'density'")
  }

  ###
  ### Prepare cluster input
  ###

  if(dims == "Comp") {
    if(!(lem %in% reducedDimNames(input))){
      stop(paste0("LEM not found, LEMs available are ", paste0(reducedDimNames(input), collapse = ", ")))
    }
    dim_dat <- reducedDim(input, lem)
    tocluster <- sampleFactors(dim_dat)
  }


  if(dims == "2d"){
    if(is.null(embedding) && is.null(xy)){
      stop("Provide an LEM and associated embedding for 2d clustering, or columns of colData()")
    }
    if(!is.null(xy)){
      if(!(xy[1] %in% names(colData(input))) || !(xy[2] %in% names(colData(input)))){
        stop("Columns not found in colData()")
      }
      tocluster <- colData(input)[,xy]
    } else {
      tocluster <- embedding(input, lem, embedding)
    }
  }

  ###
  ###
  ###


  if(method == "spectral"){
    if(is.null(k)){
      stop("Please provide k argument")
    }
    spec <- kknn::specClust(tocluster, centers = k, method = 'random-walk')
    sc_clusters <- spec$cluster
  }

  if(method == "density"){
    dist_tSNE = dist(tocluster)
    clust_tSNE = densityClust::densityClust(dist_tSNE, gaussian = T) # calculate density
    comb = as.data.frame(clust_tSNE$rho*clust_tSNE$delta) # combine rho and delta values
    colnames(comb) = "gamma"
    comb = comb[order(comb$gamma, decreasing = T), ,drop=F]
    comb$index = seq(nrow(comb))
    if (is.null(k)) {
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
      cellcut = rownames(comb)[1:k]
    }
    cellidx = which(names(clust_tSNE$rho)%in%cellcut)
    sc_clusters = densityClust::findClusters(clust_tSNE, peaks = cellidx)
    sc_clusters <- sc_clusters$clusters
  }
  colData(input)[name] <- paste0("Cluster_", sc_clusters)

  # This section is for prediction.strength() from fpc
  # It outputs the correct CBI format, however issues trying to get them to work together...
  # if(CBI){
  #   cbi_res <- vector(mode = "list", length = 5)
  #   names(cbi_res) <- c("result", "nc", "clusterlist", "partition", "clustermethod")
  #   cbi_res$clustermethod <- method
  #   cbi_res$result <- tocluster
  #   cbi_res$nc <- k
  #   cbi_res$partition <- sc_clusters
  #   clusterlist <- vector(mode = "list", length = k)
  #   clusters_cbi <- sort(unique(sc_clusters))
  #   for (i in 1:length(clusters_cbi)) {
  #     int_cluster <- clusters_cbi[i]
  #     log_vec <- sc_clusters == int_cluster
  #     names(log_vec) <- rownames(tocluster)
  #     clusterlist[[i]] <- log_vec
  #   }
  #   cbi_res$clusterlist <- clusterlist
  #   return(cbi_res)
  #
  # }

  if(cluster_stats){
    cstat <- fpc::cluster.stats(dist(tocluster), sc_clusters)
    cstat <- tibble(cstat$average.between,
                    cstat$average.within,
                    cstat$within.cluster.ss,
                    cstat$avg.silwidth,
                    cstat$pearsongamma,
                    cstat$dunn,
                    cstat$dunn2,
                    cstat$entropy,
                    cstat$wb.ratio,
                    cstat$ch)

    colnames(cstat) <- gsub("cstat", "", colnames(cstat))
    colnames(cstat) <- gsub("[[:punct:]]", "", colnames(cstat))

    # Terrible pesky argument grabbing. I do not know why this cannot be simpler...
    # tmp is a weird format, simple conversions do not work
    tmp <- mget(names(formals()),sys.frame(sys.nframe()))
    dat <- c()
    for (i in 1:length(tmp)) {
      args <- names(tmp)[i]
      val <- as.vector(unlist(tmp[args]))
      if(is.null(val)){val <- NA}
      dat <- c(dat, args, val)
    }
    dat <- (unlist(dat[3:length(dat)]))
    dat <- matrix(dat, ncol = length(dat)/2)
    colnames(dat) <- dat[1,]
    dat <- data.frame(dat)
    dat <-dat[2,,drop = F]
    cstat <- tibble(cbind(data.frame(cstat), dat))
    if(is.null(int_metadata(input)$cluster_metrics)){
      int_metadata(input)$cluster_metrics <- cstat
    } else {
      int_metadata(input)$cluster_metrics <- rbind(int_metadata(input)$cluster_metrics,cstat)

    }

  }
  return(input)
}



