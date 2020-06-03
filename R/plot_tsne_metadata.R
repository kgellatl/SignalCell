#' tSNE Plot on metadata
#'
#' This will plot information onto a 2d tsne plot
#'
#' @param input the input ex_sc.
#' @param color_by What to color points by, either "UMI_sum", or pData categorial variable, ignored if gene is provided
#' @param facet_by What to break the plots by
#' @param coords A vector of character lem, x, and u. See set_xy function.
#' #' @param title The title
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @param theme the ggplot theme
#' @param size The size of the points
#' @param ncol How many columns if faceting
#' @param legend_dot_size Size of dot in legend
#' @param text_sizes the text sizes on the plot
#' @param whether or not to shuffle the matrix to prevent last plotted factor blocking rest
#' @param xlab character x label
#' @param ylab character y label
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)

plot_tsne_metadata <- function(input,
                               color_by,
                               facet_by = NULL,
                               coords = NULL,
                               title = NULL,
                               colors = NULL,
                               # theme = "classic",
                               size = 1.5,
                               ncol = 2,
                               legend_dot_size = 1.5,
                               text_sizes = c(20,10,5,10,5,5),
                               shuffle = T,
                               xlab = "x",
                               ylab = "y"){

  if(!is.null(coords)){
    if(length(coords) == 1){
      input <- set_xy(input, lem = coords)
    } else if (length(coords) == 3) {
      input <- set_xy(input, lem = coords[1], x = as.numeric(coords[2]), as.numeric(coords[3]))
    } else {
    stop("See set_xy for how to set coords")
    }
  }

  g_Dat <- colData(input)
  g_Dat <- as_tibble(g_Dat)
  if(shuffle){
    g_Dat <- g_Dat[sample(nrow(g_Dat)),]
  } else {
    g_Dat <- g_Dat[order(g_Dat[,color_by]),]
  }
  background <- g_Dat[,c("x", "y")]

  g <- ggplot(g_Dat)
  x <- "x"
  y <- "y"
  g <- g + geom_point(data = background, aes_string(x = x, y = y), shape = 20, size = size, col = "gray")
  g <- g + geom_point(aes_string(x = x, y = y, col = color_by), shape = 20, size = size)
  if(!is.null(facet_by)){
    g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
  }
  g <- g + theme_classic()
  g <- g + theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  if(class(colData(input)[,color_by]) == "double" || class(colData(input)[,color_by]) == "integer" || class(colData(input)[,color_by]) == "numeric" ){
    g <- g +  scale_color_gradientn(colours=c('blue', 'red', 'yellow'))
  } else {
    if(all(is.na(colors)) == FALSE){
      g <- g + scale_color_manual(values = c(colors))
    }
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size=legend_dot_size)))
  if(is.null(title)){
    g <- g + labs(title = color_by)
  } else {
    g <- g + labs(title = title)
  }
  g <- g + labs(x = xlab, y = ylab)
  plot(g)
}

