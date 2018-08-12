#' @title Plot multiple combinations of sample gene expressions in a pairs plot
#'
#' @description Plots the gene expression from multiple samples in a \code{mmt} object in a grid plot with all pairs of variables.
#'
#' @param mmt (\emph{required}) A \code{mmt} list loaded with \code{\link{mt_load}}.
#' @param samples A vector of 3 or more sample names in \code{mmt} to plot on each axis. If NULL, the default, then all samples will be plotted. (\emph{Default: } \code{NULL})
#' @param label_by replace the SampleIDs plotted in the diagonal with one(!) column in the metadata.
#' @param textsize The text size of the axis titles.
#' @param pointsize The size of points in the plot(s).
#' @param linesize The size of lines in the plot(s).
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @importFrom magrittr %>%
#' @importFrom cowplot plot_grid
#' @importFrom dplyr transmute_ mutate_all
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' data("example_mmt")
#'
#' # Plot all samples (might take some time)
#' mt_plotpairs(example_mmt)
#'
#' # Plot replicates and relabel.
#' mt_plotpairs(example_mmt,samples = c("HQ180323_13","HQ180323_14"),label_by = "Replicate")
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

mt_plotpairs <- function(mmt,
                         samples   = NULL,
                         label_by  = NULL,
                         textsize  = 5,
                         linesize  = 0.5,
                         pointsize = 2){

  if(is.null(samples)){
    samples <- mmt$mtmeta[[1]]
  }
  if(!all(samples %in% mmt$mtmeta[[1]])){
    mis <- samples %w/o% mmt$mtmeta[[1]]
    stop(paste("Following sample(s) is not in the data:",paste(mis,collapse = ", ")),call. = FALSE)
  }

  ## Label by metadata instead.
  if(!is.null(label_by)){
    if (label_by %in% colnames(mmt$mtmeta)){
      label_by <- mmt$mtmeta[[label_by]]
    } else {
      stop("Your 'label_by' is not in metadata.",call. = FALSE)
    }
  } else {
    label_by <- samples
  }

  ## Prepare the data.
  dat <- mmt$mtdata[,-1] %>% {log2(. + 1)}
  rng <- range(dat)

  ## Make a blank plot
  emp <- data.frame(x = 0, y = 0)

  pblank <- ggplot(emp, aes(x,y)) +
    geom_blank() +
    theme_bw() +
    theme(plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())

  ## Iterate through the different combinations of plots.
  temp <- list()
  for (i in 1:length(samples)){
    for (j in 1:length(samples)){
      # Subset data and calc correlation (0's omitted).
      dat_sub <- dat %>%
        transmute_(x = samples[j],y = samples[i]) %>%
        mutate_all(funs(round(.,2)))

      if (i < j){
        # Remove duplicated datapoints before plotting.
        dat_sub_uni <- dat_sub %>%
        {.[!duplicated(.),]}

        # make plot
        p <- ggplot(dat_sub_uni,aes(x = x,
                           y = y)) +
          geom_point(alpha = 0.1,size = pointsize,shape = 19) +
          coord_cartesian(xlim = rng,ylim = rng) +
          geom_smooth(
            data    = dat_sub,
            method  = "gam",
            se      = F,
            formula = y ~ s(x, bs = "cs"),
            colour  = "red",
            size    = linesize) +
          geom_abline(slope = 1,intercept = 0,colour = "red",linetype = "dotted",size = linesize) +
          theme(plot.margin = margin(3,3,0,0, unit = "pt"),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA, size = 0.5, linetype = 3, color = "gray75"))
      }
      if (i == j){
        p <- pblank + geom_text(label = label_by[i], size = textsize,fontface = "bold")
      }
      if(i > j){
        # Calculate correlation (omit cases with zeros).
        r <- dat_sub %>%
        {.[apply(.,1,function(x){!any(x == 0)}),]} %>%
        {cor(.)[1,2]}
        txt <- format(c(r, 0.123456789), digits = 2)[1]
        p <- pblank + geom_text(label = txt,size = 5*r^2)
      }

      plotnr <- paste0("x_",samples[j],"y_",samples[i])
      temp[plotnr] <- list(p)
    }
  }
  ncol <- temp %>%
    length(.) %>%
    sqrt(.) %>%
    floor(.) %>%
    as.integer(.)
  cowplot::plot_grid(plotlist = temp,ncol = ncol,align = "hv")
}
