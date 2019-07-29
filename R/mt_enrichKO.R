#' mt_enrichKO
#'
#' @title Perform enrichment analysis using KEGG annotation.
#'
#' @usage mt_enrich(mmt,GeneIDs,KOs,type)
#'
#' @param mmt (\emph{required}) An object of class \code{mmt}.
#' @param GeneIDs (\emph{required}) Character vector with gene IDs of interest.
#' @param KOs (\emph{required}) String pointing to the column in mtgene containing KO identifiers.
#' @param type (\emph{required}) String refering to which KEGG annotation scheme to use. Choose from the following:
#' \itemize{
#'   \item \strong{pathway} - TBD.
#'   \item \strong{module} - TBD.
#'   \item \strong{reaction_class} - TBD.
#'   \item \strong{reaction} - TBD.
#' }
#' @param alternative (\emph{optional}) String specifying the direction of the test. Choose from:
#' \itemize{
#'   \item \strong{greater} - test for enrichment in genes of interest.
#'   \item \strong{less} - test for depletion in genes of interest.
#'   \item \strong{two.sided} - perform a two-sided test for both in genes of interest.
#' }
#'
#' @return A data.frame
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
#' # TBD.
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}
mt_enrichKO <- function(mmt,GeneIDs,type,KOs,alternative = "greater"){

  # Use all KOs identified as background universe.
  uni_genes <- mmt$mtgene[[KOs]] %>%
    ifelse(is.na(.),"",.) %>%
    lapply(function(x){
      if(x == "") x <- "" else x <- strsplit(x,split = ";")[[1]]
      return(x)
    }) %>% unlist()

  # Get the KOs from query genes.
  query_genes <- mmt$mtgene[GeneID %in% GeneIDs][[KOs]] %>%
    ifelse(is.na(.),"",.) %>%
    lapply(function(x){
      if(x == "") x <- "" else x <- strsplit(x,split = ";")[[1]]
      return(x)
    }) %>% unlist()

  # Test on type.
  grps <- lapply(KEGG_DB[[type]],"[[","KOs")

  # Perform the tests and tidy data.
  out <- t(sapply(grps,test_pathway,q_genes = query_genes,u_genes = uni_genes,alternative = alternative))
  out <- as.data.frame(out) %>% rownames_to_column(type) %>% mutate(padj = p.adjust(p,method = "bonferroni"))

  # Merge with description.
  out <- sapply(KEGG_DB[[type]],"[[","Description") %>%
  {data.frame(Description = .)} %>%
    rownames_to_column(type) %>%
    left_join(out,.,by = type) %>%
    dplyr::select(one_of(type),Description,everything())

  return(out)
}
