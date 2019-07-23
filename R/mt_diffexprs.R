#' mt_diffexprs
#'
#' @title Calculate differential expression statistics by a grouping variable.
#'
#' @usage mt_diffexprs(obj,group)
#'
#' @param obj     (\emph{required}) An object of class 'mmt' (see \code{\link{mt_load}}) or 'ampvis2' (see \code{amp_load}).
#' @param group   (\emph{required}) Character string defining the grouping variable to use for differential expression.
#' @param Genome (\emph{optional}) If provided, character string defining the grouping variable in \code{mt$mtgene} to perform genome-specific normalization (See details).
#' @param intercept  (\emph{optional}) Specify the level of 'group' to be used as reference/intercept. Default is by alphabetical order.
#' @param row_labels (\emph{optional}) Specify columns of rowdata to print on heatmap.
#' @param row_show (\emph{optional}) The number of rows to show.
#' @param order_var_by (\emph{optional}) How to order rows of heatmap. You can specify 1) "variance" or 2) one of the levels in \code{group}, which will order according to the absolute value of that level. (\emph{default}: "variance")
#'
#' @return A list with 2 elements:
#'
#' \itemize{
#'   \item \strong{heatmap} - A heatmap showing the log fold-change of each level in \code{group} relative to the intercept/baseline level specified in \code{intercept}. The \code{intercept} is therefore not possible to plot.
#'   \item \strong{DESeq} - A DESeqDataSet containing the results of the differential expression analysis. To extract results and do comparisons between groups use the \code{\link{mt_dumpDE}} function.
#' }
#'
#' @import ggplot2
#' @import ampvis2
#' @importFrom magrittr %>%
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom dplyr select rename right_join mutate
#' @importFrom tidyr gather
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactorsForMatrix DESeq results
#' @importFrom rlang syms
#'
#' @export
#'
#' @details Typically one would subset low-expression genes (use \code{\link{mt_subset}})
#' before performing the analysis, due to a low signal/noise ratio. This also drastically
#' reduces computation times and improves stability of the parameter estimation. If the
#' parameter \code{Genome} is set, a genome-specific normalization of the count matrix
#' is performed prior to standard DESeq2 analysis.
#' For more details see \href{https://peerj.com/articles/3859/}{Klingenberg & Meinicke (2017)}.
#'
#' @examples
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
#' mt_res <- mt_dumpDE(DE_mt,nom = "ELECTRODE",denom = "SUSPENTION")
#'
#' # Show the output.
#' mt_res$MAplot
#' mt_res$BOXplot
#' head(mt_res$Table)
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

mt_diffexprs <- function(obj,group,
                         Genome     = NULL,
                         intercept  = NULL,
                         row_labels = NULL,
                         order_var_by = "variance",
                         row_show   = 10){

  if (length(attributes(obj)$normalised) > 0 | !is.null(attributes(obj)$normalised)){
    stop("The DESeq2 workflow is for raw count data only and will adjust sample-effects by itself. Consider instead to remove low-count genes without normalising data to increase speed and sensitivity.")
  }
  if (class(obj) == "mmt"){

    # Ensure correct format.
    meta    <- data.frame(obj$mtmeta,stringsAsFactors = F)
    data    <- data.frame(obj$mtdata,stringsAsFactors = F)
    varmeta <- data.frame(obj$mtgene,stringsAsFactors = F)

    data    <- `rownames<-`(data,data$GeneID) %>% select(-GeneID)
  } else if (class(obj) == "ampvis2"){

    # Set some standards.
    if(is.null(tax_empty))       tax_empty       <- "best"
    if(is.null(var_label_by))    var_label_by    <- "OTU"

    # Ensure correct format.
    meta    <- obj$metadata
    data    <- obj$abund
    varmeta <- ampvis2:::amp_rename(obj,tax_empty = tax_empty)$tax %>% select(OTU,everything())

  } else {
    stop("You can only provide an object of class 'ampvis2' or 'mmt'")
  }

  ##### Input checking #####
  if(!(group %in% colnames(meta)))
    stop(paste0("group '",group,"' is not in metadata."))

  if(!is.null(row_labels)){
    if(!(any(row_labels %in% colnames(varmeta)))){
      stop(paste0("Some/all of 'row_labels' is not in the row-data."))
    } else {
      row_labels <- c(colnames(varmeta)[1],row_labels %w/o% colnames(varmeta)[1])
    }
  } else {
    row_labels <- colnames(varmeta)[1]
  }

  if(is.null(intercept)){
    levs <- as.character(sort(unique(meta[[group]])))
  } else {
    if(!(intercept %in% meta[[group]])){
      stop(paste0("intercept '",intercept,"' is not a level in '",group,"'"))
    } else {
      levs <- as.character(sort(unique(meta[[group]]))) %>% {c(intercept,(. %w/o% intercept))}
    }
  }
  meta[[group]] <- factor(meta[[group]],levels = levs)

  ##### Run the DESeq2 workflow. #####
  if (is.null(Genome)){
    dds <- DESeqDataSetFromMatrix(
      countData = data,
      colData   = meta,
      design    = as.formula(paste0("~",group)))
    mcols(dds) <- varmeta
    dds <- DESeq(dds)
  } else {
    if (!(Genome %in% colnames(varmeta))){
      stop(paste0("'",Genome,"' is not a valid column in mtgene data."))
    } else {

      # Normalize by genome.
      data_norm <- lapply(unique(varmeta$Genome),function(x){
        # Pick out genes of interest.
        wh_genes <- filter(varmeta,Genome == x)[[1]]

        # Create the normalized counts.
        tryCatch(expr = {
          as.matrix(data[wh_genes,]) %>%
            {t(t(.)/estimateSizeFactorsForMatrix(.))} %>%
            round() %>%
            `storage.mode<-`("integer")
        },error = function(e){
          stop(paste0("The estimateSizeFactorsForMatrix() fails for genome '",x,"', likely due to insufficient data."))
        })
      }) %>%
        do.call(rbind,.)

      # Reorder rows to match obj.
      data_norm <- data_norm[rownames(data),]

      # Perform DE analysis on normalized data.
      dds <- DESeqDataSetFromMatrix(
        countData = data_norm,
        colData   = meta,
        design    = as.formula(paste0("~",group)))
      mcols(dds) <- varmeta
      dds <- DESeq(dds)
    }
  }

  ##### Build heatmap #####
  # Extract results
  res <- sapply((levs %w/o% levs[1]),function(i){
    results(dds,contrast = c(group,levs[1],i))
  })

  # Calculate logFC for each level in group.
  logFC <- lapply(names(res),function(i){
    x <- res[[i]] %>%
      data.frame() %>%
      select(log2FoldChange) %>%
      {eval(parse(text = paste0("dplyr::rename(.,'",i,"' = log2FoldChange)")))} %>%
      rownames_to_column("ID")
  }) %>% {base::Reduce(f = function(x,y){merge.data.frame(x,y,by = 1)},x = .)} %>%
  {if (ncol(.) < 2)
     rownames_to_column(.,colnames(varmeta)[1])
   else
   {eval(parse(text = paste0("dplyr::rename(.,'",colnames(varmeta)[1],"' = ID)")))}}

  # Order by 'order_var_by'
  if (order_var_by == "variance"){
    if (ncol(logFC) < 3) {
      ord <- order(abs(logFC[,2]),decreasing = T)
      message("'order_var_by': not possible to order by variance with < 3 samples. Will order by absolute value instead.")
    } else {
      ord <- order(apply(logFC[,2:ncol(logFC)],1,var),decreasing = T)
    }
  } else if (order_var_by %in% meta[[group]]){
    ord <- order(abs(logFC[,order_var_by]),decreasing = T)
  } else {
    stop("'order_var_by': wrong input. Read the docs.")
  }

  # Select top 'row_show'.
  topFC <- logFC[ord,] %>%
  {if (nrow(.) <= row_show) .[1:nrow(.),] else .[1:row_show,]}

  # Merge with varmeta.
  logFC2 <- suppressWarnings(right_join(
    {varmeta %>% select(row_labels)},
    topFC,by = colnames(varmeta)[1])) %>%
  mutate(lab = paste(!!! rlang::syms(row_labels), sep = "; ")) %>%
  mutate(lab = factor(lab,levels = rev(lab))) %>%
  select(-one_of(row_labels))

  # Convert to long.
  logFC_long <- gather(logFC2,key = group,value = logFC,-lab)

  # Build plot.
  p <- ggplot(logFC_long,aes(x = group,y = lab)) +
    geom_tile(aes(fill = logFC)) +
    geom_text(aes(label = round(logFC, 1)),size = 3) +
    scale_fill_gradient2(
      low  = "blue",
      mid  = "white",
      high = "red",
      midpoint = 0,
      na.value = "grey50",
      guide = "colourbar") +
    guides(fill = guide_legend(title = paste0("logFC\n(ref: ",levs[1],")"))) +
    theme(axis.title       = element_blank(),
          panel.background = element_blank())

  ##### Return output #####
  out <- list(heatmap = p,DESeq = dds,obj = obj)
  class(out) <- "mmtDE"
  return(out)
}
