#' mt_dumpDE
#'
#' @title Extract results from DE analysis done using \code{mt_diffexprs}.
#'
#' @usage mt_dumpDE(mmtDE,num,denom)
#'
#' @param mmtDE (\emph{required}) An object of class \code{mmtDE}. See \code{\link{mt_diffexprs}}.
#' @param num (\emph{required}) Level that should be in the numerator when doing the DE analysis.
#' @param denom  (\emph{required}) Level that should be in the denominator (i.e. the reference) when doing the DE analysis. If set to "rest", the DE analysis is done as a one-vs-all test. See XXX for details.
#' @param sub_genes The genes to include in the MA and boxplot. \emph{default}: all genes)
#' @param signif The threshold for adjusted p-value to be considered significant. (\emph{default}: 0.01)
#' @param text_MAplot_size Size of text on the MA plot. (\emph{default}: 3)
#' @param ngenes_show The number of genes to show in the boxplot. (\emph{default}: 10)
#' @param order_by Order the boxplot and table by a statistic from the DE analysis. One of: (\emph{default}: "logFC")
#' @param annot_size The size of text in the MA plot. (\emph{default}: 5)
#' \itemize{
#'   \item logFC - the absolute log2 fold-change betweem \code{num} and \code{denom}.
#'   \item padj - the adjusted p-value of the DE analysis.
#' }
#' @param table_full (\emph{logical}) Output entire results table without filtering by significance. (\emph{default}: FALSE)
#' @param ... Additional arguments parsed to the \code{\link{vis_boxplot}} function.
#'
#' @return A list with 3 elements:
#' \itemize{
#'   \item \strong{MAplot} - A MAplot (ggplot-object) with the gene expression plotted against the log2 fold-change.
#'   \item \strong{BOXplot} - A faceted boxplot showing the expression of the most interesting genes, grouped according to the DE analysis design.
#'   \item \strong{Table} - A table with the variable metadata and associated statistics from the DE analysis.
#' }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom DESeq2 results
#' @importFrom dplyr mutate arrange right_join filter select
#'
#' @export
#'
#' @details
#'
#' @examples
#'
#' \dontrun{
#' # Load data.
#' data("example_mmt")
#'
#' # Compute the statistics.
#' DE_mt <- mt_diffexprs(example_mmt,
#'   group      = "Type",
#'   row_labels = c("GeneID","product"),
#'   intercept  = "ANAMMOX")
#'
#' # Extract results for given levels.
#' mt_res <- mt_dumpDE(DE_mt,num = "ELECTRODE",denom = "SUSPENTION")
#'
#' # Show the output.
#' mt_res$MAplot
#' mt_res$BOXplot
#' head(mt_res$Table)
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}

mt_dumpDE <- function(mmtDE,num,denom,
                      sub_genes   = NULL,
                      signif      = 0.01,
                      text_MAplot_size  = 3,
                      ngenes_show = 10,
                      order_by    = "logFC",
                      table_full  = F,
                      annot_size  = 5,
                      ...){

  if(class(mmtDE) != "mmtDE"){
    stop("You need to provide data of class 'mmtDE' as generated using the mt_diffexprs function.",call. = F)
  }

  # Get group and reference level.
  group <- mmtDE$DESeq@design %>% {as.character(.)[2]}
  levs  <- {mmtDE$DESeq[[group]]} %>% levels()

  if (denom != "rest"){
    if (!all(c(num,denom) %in% mmtDE$DESeq[[group]])){
      stop(paste0("'num' and 'denom' must be in levels of the grouping variable '",group,"'"),call. = F)
    }
  }

  # Setup contrasts correctly.
  if (denom == "rest"){
    if (num == levs[1]){
      ctrst <- c(0,rep(-1/(length(levs)-1),length(levs)-1))
    } else {
      ctrst <- c(0,rep(-1/(length(levs)-2),length(levs)-1))
      ctrst[which(levs == num)] <- 1
    }
    denom <- paste0(levs[-1],collapse = "+")
  } else {
    ctrst <- c(group,num,denom)
  }

  # Do the analysis.
  res <- results(mmtDE$DESeq,contrast = ctrst)
  res <- res %>%
  {cbind.data.frame(ID = rownames(.),.)}
  if(!is.null(sub_genes)) res <- filter(res,ID %in% sub_genes)
  res <- res %>%
  mutate(Significant = factor(ifelse(padj < signif, "Yes", "No"),levels = c("Yes","No")),
         ER          = factor(ifelse(log2FoldChange > 0,num,denom),levels = c(num,denom)))

  ### MA PLOT ##################################################################
  # Count the number of differentially expressed.
  labs <- subset(res,select = c(Significant,ER)) %>%
    table %>% .[1,c(num,denom)] %>% {paste(names(.),.,sep = "\nn = ")}

  # plot.
  MAplot <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
    geom_hline(yintercept = 0, color = "darkred", lty = 2) +
    geom_point(size = 1) +
    scale_x_log10() +
    scale_color_manual(na.value = "black", values = c("red", "black")) +
    theme_classic() +
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          legend.position = "none") +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.2, label = labs[1], size = annot_size) +
    annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -.2, label = labs[2], size = annot_size) +
    ylab("Fold Change (log2)") +
    xlab("Expression (base mean)")

  ### Boxplot ##################################################################
  # Select top ngenes_show for plotting.
  if (order_by %in% c("logFC","padj")){
    res <- res %>% {
      switch(order_by,
             logFC = arrange(.,desc(abs(log2FoldChange))),
             padj  = arrange(.,padj))
    }
  } else {
    stop("order_by: please provide correct argument.",call. = FALSE)
  }

  # Build plot.
  wh <- res$ID[1:ngenes_show] %>% as.character()
  p <- vis_boxplot(obj      = mmtDE$obj,
                   group_by = as.character(mmtDE$DESeq@design)[2],
                   sort_by  = "row_show",
                   row_show = wh,
                   ...)

  ctrst_name <- paste(paste("numer:",num),paste("denom:",denom),sep = "\n")
  BOXplot <- p +
    geom_point(aes(x=Inf, y=Inf,shape = ctrst_name),alpha = 0) +
    guides(shape = guide_legend(title = "Contrast"))

  ### Table ####################################################################
  if(class(mmtDE$obj) == "mmt"){
    out <- suppressWarnings(
      right_join(mmtDE$obj$mtgene,res,by = c("GeneID" = "ID")))
  } else if (class(mmtDE$obj) == "ampvis2"){
    out <- suppressWarnings(
      right_join(mmtDE$obj$tax,res,by = c("OTU" = "ID")))
  }

  if (!table_full){
    out <- filter(out,Significant == "Yes")
  }

  out <- select(out,-Significant,-ER)

  ### Build output #############################################################
  return(list(MAplot = MAplot,BOXplot = BOXplot,Table = out))
}
