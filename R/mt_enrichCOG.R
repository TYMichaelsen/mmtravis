#' mt_enrichCOG
#'
#' @title Perform enrichment analysis using COG categories.
#'
#' @usage mt_enrich(mmt,GeneIDs,COGs)
#'
#' @param mmt (\emph{required}) An object of class \code{mmt}.
#' @param GeneIDs (\emph{required}) Character vector with gene IDs of interest.
#' @param COGs (\emph{required}) String pointing to the column in mtgene containing COG identifiers.
#' @param alternative (\emph{optional}) String specifying the direction of the test. Choose from:
#' \itemize{
#'   \item \strong{greater} - test for enrichment in genes of interest.
#'   \item \strong{less} - test for depletion in genes of interest.
#'   \item \strong{two.sided} - perform a two-sided test for both in genes of interest.
#' }
#' @param show_p (\emph{optional}) Numeric indicating the threshold p-value for plotting. Default: 1
#' @param size_p (\emph{optional}) Numeric indicating the size of plotted p-values. Default: 1
#' @param unannotated (\emph{optional}) Should the data for unannotated genes be included in plot? Default: \code{TRUE}.
#' @return A list with following:
#' \itemize{
#'   \item \strong{table} - A table containing the results from the enrichment analysis.
#'   \item \strong{plot} - A ggplot barplot showing reduction/enrichment of each entry in \code{type} in percent, relative to the percent found in the entire dataset.
#' }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate left_join select one_of everything filter starts_with group_by arrange mutate_at distinct ungroup
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather spread
#' @importFrom rlang sym
#'
#' @export
#'
#' @details
#'
#' @examples
#'
#' \dontrun{
#' # TBD.
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}
mt_enrichCOG <- function(mmt,GeneIDs,COGs,alternative = "greater",show_p = 1,size_p = 1,unannotated = T){
  type = "Category"

  # Use all COGs identified as background universe.
  uni_genes <- mmt$mtgene[[COGs]] %>%
    ifelse(is.na(.),"",.) %>%
    ifelse(. == "","None",.)

  # Get the COGs from query genes.
  query_genes <- mmt$mtgene[GeneID %in% GeneIDs][[COGs]] %>%
    ifelse(is.na(.),"",.) %>%
    ifelse(. == "","None",.)

  # Test on type. Remove the entires in type that are missing from query genes.
  grps  <- lapply(COG_DB[[type]],"[[","COGs")
  keeps <- lapply(grps,`%in%`,query_genes) %>% sapply(any)
  grps  <- grps[keeps]

  # Perform the tests and tidy data.
  tabout <- t(sapply(grps,test_pathway,q_genes = query_genes,u_genes = uni_genes,alternative = alternative))
  tabout <- as.data.frame(tabout) %>% rownames_to_column(type) %>% mutate(padj = p.adjust(p,method = "bonferroni"))

  # Merge with description.
  tabout <- sapply(COG_DB[[type]],"[[","Description") %>%
  {data.frame(Description = .)} %>%
    rownames_to_column(type) %>%
    left_join(tabout,.,by = type) %>%
    dplyr::select(one_of(type),Description,everything())

  # Make the plot.--------------------------------------------------------------
  enrich_pct <- filter(tabout,INuniverse > 0 & padj <= show_p)
  if(!unannotated){
    enrich_pct <- filter(enrich_pct,!!sym(type) != "None")
  }
  if (nrow(enrich_pct) < 1){
    p_enrich <- paste0("There is no data to plot!")
  } else {
    enrich_pct <- enrich_pct %>%
      mutate(base   = INuniverse/Nuniverse*100) %>%
      mutate(differ = (INquery/Nquery)*100 - base) %>%
      mutate(pct_uni     = ifelse(differ < 0,base + differ,base),
             pct_overlap = ifelse(differ < 0,-differ,0),
             pct_added   = ifelse(differ < 0,0,differ)) %>%
      select(!!sym(type),Description,differ,starts_with("pct_")) %>%
      gather(key = "variable",value = "value",-!!sym(type),-Description,-differ) %>%
      group_by(!!sym(type)) %>%
      mutate(total   = sum(value)) %>%
      arrange(differ) %>%
      ungroup() %>%
      mutate_at(.vars = type,.funs = list(~factor(.,levels = unique(.))))

    # Fetch p-values.
    pvals <- tabout %>%
      mutate(padj = ifelse(padj > 1,"",ifelse(padj < 0.001,"<0.001",round(padj,3)))) %>%
      mutate_at(.vars = type,.funs = list(~factor(.,levels = levels(enrich_pct[[type]])))) %>%
      left_join(select(enrich_pct,!!sym(type),"total"),.,by = type) %>%
      distinct()

    # Position of "|".
    pos_point <- enrich_pct %>%
      spread(key = "variable",value = "value") %>%
      mutate(PosX = pct_uni+pct_overlap) %>%
      mutate(PosX = ifelse(PosX < max(enrich_pct$total)*0.01/2,max(enrich_pct$total)*0.01/2,PosX))

    cols <- c("pct_added" = "#F8766D", "pct_overlap" = "#00BFC4", "pct_uni" = "lightgrey")

    # Plot.
    p_enrich <- ggplot(enrich_pct,aes_string(x = type,y = "value",fill = "variable")) +
      geom_bar(stat = "identity",width = .8) +
      coord_flip(clip = "off") +
      geom_text(aes_string(label = "padj",x = type,y = "total + max(total)*.01"),data = pvals,
                inherit.aes = F,hjust = 0,vjust = .3,size = size_p) +
      #geom_point(aes(x = x,y = y,shape = "All genes"),data = data.frame(x = Inf,y = -Inf),inherit.aes = F) +
      geom_tile(aes_string(x = type,y = "PosX"),data = pos_point,width = 1,height = max(enrich_pct$total)*0.005,inherit.aes = F) +
      scale_fill_manual(
        values = cols,
        breaks = c("pct_added","pct_overlap"),
        labels = c("Enriched", "Reduced")) +
      #scale_shape_manual(
      #  values = c("All genes" = "|")) +
      scale_y_continuous(
        expand = expand_scale(mult = c(0,.1)),
        limits = c(0,max(enrich_pct$total)*1.06)) +
      labs(
        y = "Percentage of genes [%]",
        x = "",
        fill = "") +
      theme_classic() +
      theme(legend.spacing.y = unit(0, "pt"))
  }
  return(list(table = tabout,plot = p_enrich))
}
