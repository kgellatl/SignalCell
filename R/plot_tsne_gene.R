#' tSNE Plot on a gene or gene
#'
#' This will plot gene information onto a 2d tsne plot
#'
#' @param input The input data
#' @param gene the gene or genes to plot
#' @param facet_by If ONE gene is being plotted it can be faceted. You cannot facet multiple gene
#' @param assay The assay to operate on. Will default to get_def_assay(input)
#' @param coords A vector of character lem, x, and u. See set_xy function.
#' @param title The title
#' @param colors Colors for the cells
#' @param theme the ggplot theme
#' @param size The size of the points
#' @param ncol controls the number of columns if faceting
#' @param text_sizes the text sizes on the plot
#' @param xlab character x label
#' @param ylab character y label
#' @export
#' @details
#' Utilize information stored in colData to control the plot display.
#' @examples
#' plot_tsne_gene(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)
#'
#'
plot_tsne_gene <- function(input,
                           gene,
                           facet_by = NULL,
                           assay = NULL,
                           coords = NULL,
                           title = NULL,
                           colors = c("gray", 'blue', 'red', 'yellow'),
                           # theme = "classic",
                           size = 1.5,
                           ncol = 2,
                           text_sizes = c(20,10,5,10,5,5),
                           xlab = "x",
                           ylab = "y") {

  if(length(gene) > 1 && !(is.null(facet_by))){
    stop("Cannot facet multiple genes.")
  }

  if(!(gene %in% rownames(input))){
    stop("Gene not found")
  }

  if(is.null(assay)){
    def_assay <- get_def_assay(input)
    input_mat <- assay(input, def_assay)[gene,]
    assay_name <- def_assay
  } else {
    if(!(assay%in%names(input@assays))){
      stop(paste0("Assay not found, assays available are, ", names(input@assays)))
    }
    input_mat <- assay(input, assay)[gene,]
    assay_name <- assay
  }

  if(!is.null(coords)){
    if(length(coords) == 1){
      input <- set_xy(input, lem = coords)
    } else if (length(coords) == 3) {
      input <- set_xy(input, lem = coords[1], x = as.numeric(coords[2]), as.numeric(coords[3]))
    } else {
      stop("See set_xy for how to set coords")
    }
  }

  input_mat <- as.matrix(input_mat)
  input_mat <- as.data.frame(input_mat)
  if(length(gene) > 1){
    input_mat <- t(input_mat)
  }
  input_mat <- cbind(input_mat, colData(input))
  colnames(input_mat)[1:length(gene)] <- gene
  input_mat <- as_tibble(input_mat)

  input_mat <- pivot_longer(input_mat, cols = 1:length(gene), names_to = "gene", values_to = "Expression")

  g_Dat <- input_mat
  background <- g_Dat[,c("x", "y")]

  g_Dat <- g_Dat[which(g_Dat$Expression > 0),]

  g <- ggplot(g_Dat)
  g <- g + geom_point(data = background, aes(x = x, y = y), shape = 20, size = size, col = "gray")
  g <- g +  scale_color_gradientn(colours=colors)
  g <- g +  geom_point(aes(x=x, y=y, col = Expression), shape=20, size = size)
  if(length(gene > 1)){
    g <- g +  facet_wrap(~gene, ncol = ncol)
  }
  if(!is.null(facet_by)){
    g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol)
  }
  g <- g +  theme_classic()
  g <- g +  theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
  g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  if(!is.null(title)){
    g <- g +  labs(title = title)
  } else {
    if(!is.null(facet_by)){
      g <- g +  labs(title = gene)
    }
  }
  g <- g + labs(col = assay_name, x = xlab, y = ylab)
  return(g)
}

