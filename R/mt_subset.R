#' Subset mt objects based on sample metadata
#'
#' Subsets the data in ampvis2 objects based on metadata and returns the subsetted object.
#'
#' @usage mt_subset(data, ...)
#'
#' @param mt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param sub_genes (\emph{optional}) A string specifying a logical row subset operation on the mtgene dataframe in the mt object parsed to \code{subset()}.
#' @param sub_samples (\emph{optional}) A string specifying a logical row subset operation on the mtmeta dataframe in the mt object parsed to \code{subset()}.
#' @param minreads Minimum number of reads pr. gene. Genes below this value will be removed. (\emph{default:} \code{0})
#'
#' @return A modifed mt object
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @export
#'
#' @details The subset is performed on the mtmeta/mtgene data by \code{subset()} and the whole object is then adjusted accordingly.
#'
#' @example Some example here.
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

mt_subset <- function(mt,sub_genes = NULL,sub_samples = NULL,minreads = 0) {

  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(mt$meta) %>% as.numeric()
  ngenesbefore   <- nrow(mt$count) %>% as.numeric()

  #remove genes below minreads.
  wh <- rowSums(mt$count[,-1]) >= minreads
  mt$count <- mt$count[wh,,drop = F]
  mt$gene  <- if(!is.null(mt$gene)) mt$gene <- mt$gene[wh,,drop = F]

  # Subset samples.
  if (!is.null(sub_samples)){
    wh <- tryCatch(expr = {
      mt$meta %>% filter(eval(parse(text = sub_samples)))
    },error = function(e){
      stop("The provided 'sub_samples' string is not meaningfull for the metadata.")
    }) %>% {mt$meta$SampleID %in% .$SampleID}
    mt$meta  <- mt$meta[wh,,drop = F]
    mt$count <- mt$count[,c(TRUE,wh),drop = F]
    if(!is.null(mt$stat)) mt$stat <- mt$stat[wh,,drop = F]
  }

  # Subset genes.
  if (!is.null(mt$gene)){
    if(!is.null(sub_genes)){
      wh <- tryCatch(expr = {
        mt$gene %>% filter(eval(parse(text = sub_genes)))
      },error = function(e){
        stop("The provided 'sub_genes' string is not meaningfull for the gene data.")
      }) %>% {mt$gene$GeneID %in% .$GeneID}
      mt$count <- mt$count[wh,,drop = F]
      mt$gene  <- mt$gene[wh,,drop = F]
    }
  } else {
    stop("There is no gene data available for this mt object.")
  }

  #remove genes below minreads.
  wh <- rowSums(mt$count[,-1]) >= minreads
  mt$count <- mt$count[wh,,drop = F]
  mt$gene  <- if(!is.null(mt$gene)) mt$gene <- mt$gene[wh,,drop = F]

  # Report the removed samples/genes.
  nsamplesafter <- nrow(mt$meta) %>% as.numeric()
  ngenesafter   <- nrow(mt$count) %>% as.numeric()

  message(paste(nsamplesbefore - nsamplesafter, "samples and",
                ngenesbefore - ngenesafter, "genes have been filtered \nBefore:",
                nsamplesbefore, "samples and", ngenesbefore, "genes\nAfter:",
                nsamplesafter, "samples and", ngenesafter, "genes"))
  return(mt)
}
