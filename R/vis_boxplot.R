#' Boxplot
#'
#' Generates boxplots of genes/OTUs.
#'
#' @usage vis_boxplot(obj,group_by)
#'
#' @param obj (\emph{required}) A data object of class 'ampvis2' (see \code{\link[ampvis2]{amp_load}}) or class 'mmt' (see \code{\link{mt_load}}).
#' @param group_by (\emph{required}) Group the samples by a variable in the metadata.
#' @param boxgroup_by Group the samples within each boxplot by a variable in the metadata.
#' @param order_group A vector to order the groups by.
#' @param order_boxgroup A vector to order the box-groups by.
#' @param row_show The number of variables to show, or a vector of gene/OTU names. (\emph{default:} \code{20})
#' @param title_newline Split the facet headers by newline. (\emph{default:} \code{FALSE})
#' @param title_size Size of text in facet headers. (\emph{default:} \code{5})
#' @param plot_flip (\emph{logical}) Flip the axes of the plot axis. (\emph{default:} \code{FALSE})
#' @param plot_log (\emph{logical}) Log10-scale the plot. (\emph{default:} \code{FALSE})
#' @param adjust_zero Keep values of 0 in calculations by adding this value. (\emph{default:} \code{NULL})
#' @param point_size The size of points. (\emph{default:} \code{1})
#' @param sort_by Sort the boxplots by \code{"median"}, \code{"mean"} or \code{"total"}. (\emph{default:} \code{"median"})
#' @param plot_type Plot type. \code{"boxplot"} or \code{"point"}. (\emph{default:} \code{"boxplot"})
#' @param normalise Can take one of following arguments: (\emph{default:} \code{"libsize"})
#' \itemize{
#'    \item \code{"total"}: Normalise the read counts to be in percent per sample.
#'    \item \code{"libsize"}: Normalise the read counts to adjust for gene dispersion and total read counts per sample. See \link[DESeq2]{estimateSizeFactorsForMatrix} for details.
#'    \item \code{"none"}: No normalisation.
#'    }
#' @param detailed_output (\emph{logical}) Return additional details or not. If \code{TRUE}, it is recommended to save to an object and then access the additional data by \code{View(object$data)}. (\emph{default:} \code{FALSE})
#' @param row_labels (\strong{mmtravis only}) Additional variables in genedata to plot (\code{mmt}).
#' @param tax_class (\strong{ampvis2 only}) Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param tax_aggregate (\strong{ampvis2 only}) The taxonomic level to aggregate the OTUs. (\emph{default:} \code{"Genus"})
#' @param tax_empty (\strong{ampvis2 only}) How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible.
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#'
#' @return A ggplot2 object. If \code{detailed_output = TRUE} a list with a ggplot2 object and additional data.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate right_join rename one_of group_by summarise arrange filter
#' @importFrom tibble rownames_to_column
#' @importFrom rlang syms
#' @importFrom tidyr gather
#' @importFrom data.table as.data.table
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#'
#' @export
#'
#' @details The plot is ordered as one would read a book, from left to right and top to bottom.
#'
#' @examples
#'
#' \dontrun{
#' data("example_mmt")
#'
#' vis_boxplot(example_mmt,
#'   group_by    = "Type",
#'   boxgroup_by = "Organism",
#'   row_show    = 10,
#'   row_labels  = "product",
#'   title_newline = T)
#' }
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

vis_boxplot <- function(obj,group_by,
                        boxgroup_by = NULL,
                        sort_by = "median",
                        normalise = "libsize",
                        plot_type = "boxplot",
                        point_size = 1,
                        row_show = 20,
                        order_group = NULL,
                        order_boxgroup = NULL,
                        title_newline = FALSE,
                        title_size = 5,
                        plot_flip = FALSE,
                        plot_log = FALSE,
                        adjust_zero = NULL,
                        detailed_output = FALSE,
                        # mmtravis-specific.
                        row_labels = NULL,
                        # Ampvis-specific.
                        tax_class = NULL,
                        tax_add   = NULL,
                        tax_empty = NULL,
                        tax_aggregate = NULL){

  ### Checking data format #####################################################
  "%w/o%" <- function(x, y) x[!x %in% y]

  if (class(obj) == "mmt"){

    ### Set/check options. ###
    if (!is.null(tax_empty)) message("'tax_empty' has no effect on 'mmt' objects, ignore it.")
    if (!is.null(tax_aggregate)) message("'tax_aggregate' has no effect on 'mmt' objects, ignore it.")
    if (!is.null(tax_add)) message("'tax_add' has no effect on 'mmt' objects, ignore it.")
    if (!is.null(tax_class)) message("'tax_add' has no effect on 'mmt' objects, ignore it.")

    id <- "GeneID"
    if(is.null(row_labels)){
      row_labels <- "GeneID"
    } else {
      row_labels <- c(id,row_labels %w/o% id)
    }

    ### Ensure correct format. ###
    meta    <- obj$mtmeta
    data    <- obj$mtdata
    varmeta <- obj$mtgene

    if(!all(row_labels %in% colnames(varmeta))) stop("row_labels: not in variable metadata.",call. = FALSE)

    data    <- `rownames<-`(data,data$GeneID) %>% select(-GeneID)

    # Check for normalization.
    if(!is.null(attributes(obj)$normalised)){
      allreadyNormalised <- T
      if (normalise != "none"){
        warning("The data has already been normalised by mt_subset. You should set normalise = 'none' to avoid losing information about the original data of which the provided data is a subset of.", call. = FALSE)
      }
    }

    ### Prepare for plotting ###
    data_wide <- suppressMessages(right_join(
      {varmeta %>% select(row_labels)},
      {data %>% rownames_to_column(id)})) %>%
      mutate(.,Display = paste(!!! rlang::syms(row_labels), sep = ifelse(title_newline,"\n","; "))) %>%
      select(-one_of(row_labels %w/o% id)) %>%
      {eval(parse(text = paste0("dplyr::rename(.,ID = '",id,"')")))} %>%
      select(ID,Display,everything())

  } else if (class(obj) == "ampvis2"){

    ### Set/check options. ###
    if (!is.null(row_labels)) message("'row_labels' has no effect on 'ampvis2' objects, use 'tax_aggregate' and 'tax_add' instead.")

    if(is.null(tax_empty))     tax_empty     <- "best"
    if(is.null(row_labels))    row_labels    <- "Genus"
    if(is.null(tax_aggregate)) tax_aggregate <- "Genus"
    id <- "OTU"

    # tax_add and tax_aggregate can't be the same
    if(!is.null(tax_aggregate) & !is.null(tax_add)) {
      if(tax_aggregate == tax_add) {
        stop("tax_aggregate and tax_add cannot be the same", call. = FALSE)
      }
    }

    ### Ensure correct format. ###
    meta    <- obj$metadata
    data    <- obj$abund
    varmeta <- ampvis2:::amp_rename(obj,
      tax_empty = tax_empty,
      tax_class = tax_class,
      tax_level = tax_aggregate)$tax %>% select(OTU,everything())

    # Check for normalization.
    if(isTRUE(attributes(obj)$normalised)){
      allreadyNormalised <- TRUE
      if (normalise != "none"){
        warning("The data has already been normalised by either amp_subset_samples or amp_subset_taxa. You should set normalise = 'none' to avoid losing information about the original data of which the provided data is a subset of.", call. = FALSE)
      }
    }

  ### Prepare for plotting (including aggregation) ###
    data_wide <- cbind.data.frame(varmeta,data) %>%
      mutate(.,Display = paste(!!! rlang::syms(row_labels), sep = ifelse(title_newline,"\n","; "))) %>%
      select(-one_of(colnames(varmeta) %w/o% tax_aggregate)) %>%
      tidyr::gather(key = Sample, value = Abundance,-Display,-one_of(tax_aggregate)) %>%
      as.data.table() %>%
      {eval(parse(text = paste0("dplyr::rename(.,ID = '",tax_aggregate,"')")))} %>%
      .[,sum(Abundance),by = list(ID,Display,Sample)] %>%
      dplyr::rename(Abundance = V1) %>%
      tidyr::spread(key = Sample,value = Abundance)

  } else {
    stop("You can only provide an object of class 'ampvis2' or 'mmt'",call. = FALSE)
  }

  ### Sanity check of options ##################################################
  if(!all(group_by %in% colnames(meta))) stop("group_by: not in metadata.",call. = FALSE)
  if(!is.null(boxgroup_by)){
    if(!all(boxgroup_by %in% colnames(meta))) stop("boxgroup_by: not in metadata.",call. = FALSE)
    if(boxgroup_by == group_by) stop("You cannot set boxgroup_by and group_by to the same.",call. = FALSE)
  }


  ### Prepare data for plotting ################################################
  # Normalise.
  if (normalise == "total" | normalise == "none"){
  data_long <- data_wide %>%
    tidyr::gather(key = Sample, value = Abundance, -Display,-ID) %>%
    as.data.table() %>%
    group_by(Sample) %>%
    {if (normalise == "none") . else mutate(.,Abundance = prop.table(Abundance)*100)}
  } else if (normalise == "libsize"){
    x <- data_wide %>%
      select(-Display,-ID) %>%
      as.matrix() %>%
      {t(t(.)/DESeq2::estimateSizeFactorsForMatrix(.))} %>%
      as.data.frame()
    data_long <- cbind.data.frame(select(data_wide,ID,Display),x) %>%
      tidyr::gather(key = Sample, value = Abundance, -Display,-ID)
  }

  ### ORDERING ###
  t1 <- sort_by == "row_show"
  if (is.numeric(row_show)){
    t2 <- F
  } else {
    if (all(row_show %in% data_wide$ID)){
      t2 <- T
    } else {
      stop("row_show: you are not pointing to valid row IDs. This should be a taxonomic label corresponding to tax_aggregate for ampvis2 and Gene IDs for mmtravis",call. = FALSE)
    }
  }

  if (t1 & t2){
    wh <- row_show
  } else if (t1 & !t2) {
    stop("sort_by: you cannot use 'row_show' in sort_by, when row_show itself is a number.",call. = FALSE)
  } else {

    TotalCounts <- group_by(data_long,ID) %>%
      summarise(Median = median(Abundance), Total = sum(Abundance), Mean = mean(Abundance))

    if        (sort_by == "median"){wh <- TotalCounts %>% arrange(desc(Median)) %>% as.data.frame() %>% .$ID
    } else if (sort_by == "mean"){wh   <- TotalCounts %>% arrange(desc(Mean)) %>% as.data.frame() %>% .$ID
    } else if (sort_by == "total"){wh  <- TotalCounts %>% arrange(desc(Total)) %>% as.data.frame() %>% .$ID
    } else stop("sort_by: please specify a valid argument.",call. = FALSE)

    if(!t1 & t2){
      # Keep only rows from "row_show", but retain ordering.
      wh <- wh[wh %in% row_show]
    } else if (!t1 & !t2){
      # retain only top row_show.
      wh <- wh[1:row_show]
    }
  }

  # Subset data.
  data_long_sort <- data_long %>%
    filter(ID %in% wh) %>%
    mutate(ID = factor(ID,levels = wh))

  # Merge with grouping variable.
  suppressWarnings(
    if (group_by != colnames(meta)[1]){
      if (length(group_by) > 1){
        grp <- data.frame(Sample = meta[[1]], Group = apply(meta[,group_by], 1, paste, collapse = " "))
      } else{
        grp <- data.frame(Sample = meta[[1]], Group = meta[,group_by])
      }
      data_long_sort$Group <- grp$Group[match(data_long_sort$Sample, grp$Sample)]
    } else{ data_long_sort <- data.frame(data_long_sort, Group = data_long_sort$Sample)}
  )

  # Merge with box-grouping variable.
  suppressWarnings(
    if (!is.null(boxgroup_by)){
      if (boxgroup_by != colnames(meta)[1]){
        if (length(boxgroup_by) > 1){
          grp <- data.frame(Sample = meta[[1]], Group = apply(meta[,boxgroup_by], 1, paste, collapse = " "))
        } else{
          grp <- data.frame(Sample = meta[[1]], Group = meta[,boxgroup_by])
        }
        data_long_sort$BoxGroup <- grp$Group[match(data_long_sort$Sample, grp$Sample)]
      }
    }
  )

  # Add a small constant to handle ggplot2 removal of 0 values in log scaled plots
  if(!is.null(adjust_zero)){
    data_long_sort$Abundance[data_long_sort$Abundance == 0] <- adjust_zero
  }

  # Order group.
  if(!is.null(order_group)){
    data_long_sort$Group <- factor(data_long_sort$Group, levels = rev(order_group))
  }

  # Order boxgroup.
  if(!is.null(order_boxgroup)){
    data_long_sort$BoxGroup <- factor(data_long_sort$BoxGroup, levels = rev(order_boxgroup))
  }

  ### PLOTTING #################################################################

  if (!is.null(data_long_sort$BoxGroup)){
    p <- ggplot(data_long_sort, aes(x = BoxGroup, y = Abundance,fill = Group)) +
      facet_wrap(~ Display) +
      theme(axis.text.x  = element_text(angle = 90,vjust = 0.5,hjust = 1))
  } else {
    p <- ggplot(data_long_sort, aes(x = Group, y = Abundance,fill = Group)) +
      facet_wrap(~ Display) +
      theme(axis.text.x  = element_blank())
  }

  p <- p +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_line(color = "grey90"),
          strip.text.x       = element_text(size = title_size)) +
    xlab("")

  switch(normalise,
         none    = if (allreadyNormalised) p <- p + ylab("Normalised Abundance") else p <- p + ylab("Raw counts [#]"),
         total   = p <- p + ylab("Relative Abundance [%]"),
         libsize = p <- p + ylab("Normalized counts [#]"))

  if (plot_flip == TRUE){ p <- p + coord_flip()}

  if (plot_type == "point"){ p <- p + geom_point(aes(colour = Group),size = point_size)}
  if (plot_type == "boxplot"){p <- p + geom_boxplot(outlier.size = point_size)}
  if (plot_log == TRUE){ p <- p + scale_y_log10()}

  p <- p +  guides(fill = guide_legend(title = paste(group_by,collapse = " - ")))
  if (detailed_output) {
    return(list(plot = p, data = data_long_sort))
  } else
    return(p)
}
