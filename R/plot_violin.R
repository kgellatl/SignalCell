#' Plot Violin
#'
#' This will plot gene expression via violin plot
#'
#' @param input The input data
#' @param gene the gene or genes to plot
#' @param color_by colData() column to facet plot by
#' @param facet_by If ONE gene is being plotted it can be faceted. You cannot facet multiple genes
#' @param assay The assay to operate on. Will default to get_def_assay(input)
#' @param title The title
#' @param jitter Whether or not to add jitter points to violin
#' @param plot_mean whether or not to plot the mean value along with a secondary axis
#' @param number_labels Add counts and fraction expression to plot
#' @param sig number of decimals for number_labels
#' @param size The size of the points
#' @param ncol controls the number of columns if faceting
#' @param text_sizes the text sizes on the plot
#' @export
#' @details
#' Utilize information stored in colData to control the plot display.
#' @examples
#' plot_violin(ex_sc_example, gene = "Tnf", title = "Tnf over Time", facet_by = "Timepoint", density = TRUE)
#'
#'
plot_violin <- function(input,
                        gene,
                        color_by,
                        facet_by = NULL,
                        assay = NULL,
                        title = NULL,
                        # theme = "classic",
                        jitter = T,
                        plot_mean = T,
                        number_labels = T,
                        sig = 2,
                        size = 1,
                        ncol = 2,
                        text_sizes = c(20,10,5,10,5,5)) {

  if(length(gene) > 1 && !(is.null(facet_by))){
    stop("Cannot facet multiple genes.")
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

  if(!is.null(facet_by)){
    sum_groups <- c(color_by, facet_by)
  } else {
    if(length(gene) >1 ){
      sum_groups <- c(color_by, "gene")
    } else {
      sum_groups <- color_by
    }
  }

  frac_finder <- function(x){
    length(which(x>0))/length(x)
  }

  summary_dat = function(x) {
    g_Dat %>%
      group_by_at(sum_groups) %>%
      summarise(mean = mean(Expression), n = n(), frac = frac_finder(Expression), .groups = "rowwise")
  }

  if(number_labels || plot_mean){
    summary_data <- summary_dat(g_Dat)
    summary_data$frac <- round(summary_data$frac, sig)
    summary_data$num_pos <- -max(g_Dat$Expression)/10
    summary_data$frac_pos <- -max(g_Dat$Expression)/5
  }

  g <- ggplot(g_Dat)
  g <- g +  geom_violin(aes_string(x = color_by, y = "Expression", col = color_by), size = size,  trim = T, scale = "width")
  if(number_labels){
    g <- g +  geom_text(data = summary_data, aes_string(x = color_by, y = "num_pos", label = "n"), size = size*2)
    g <- g +  geom_text(data = summary_data, aes_string(x = color_by, y = "frac_pos", label = "frac"), size = size*2)
  }
  if(jitter){
    g <- g + geom_jitter(aes_string(x=color_by, y="Expression", col = color_by), width = 0.2, size = size)
  }
  if(plot_mean){
    g <- g +  geom_point(data = summary_data, aes_string(x = color_by, y = "mean"), size = size*2)
  }
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
    g <- g +  labs(title= title)
  } else {
    if(!is.null(facet_by)){
      g <- g +  labs(title= gene)
    }
  }
  g <- g + labs(y = assay_name)
  return(g)
}


