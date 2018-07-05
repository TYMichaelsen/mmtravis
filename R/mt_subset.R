#' Subset mt objects based on sample metadata
#'
#' @title Normalise and subset the data in \code{mmt} objects based on variable/sample metadata.
#'
#' @usage mt_subset(data, ...)
#'
#' @param mt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param sub_genes (\emph{optional}) A string specifying a logical row subset operation on the mtgene dataframe in the mt object parsed to \code{subset()}.
#' @param sub_samples (\emph{optional}) A string specifying a logical row subset operation on the mtmeta dataframe in the mt object parsed to \code{subset()}.
#' @param minreads Minimum number of reads pr. gene. Genes below this value will be removed. (\emph{default:} \code{0})
#' @param normalise Transform the read counts AFTER reads have been removed by the minreads argument. (\emph{default:} \code{"none"})
#' \itemize{
#'    \item \code{"total"}: Normalise the read counts to be in percent per sample.
#'    \item \code{"libsize"}: Normalise the read counts to adjust for gene dispersion and total read counts per sample. See \link[DESeq2]{estimateSizeFactorsForMatrix} for details.
#'    \item \code{"none"}: No normalisation.
#'    }
#'
#' @return A modifed mt object
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @export
#'
#' @details The subset is performed on the mtmeta/mtgene data by \code{subset()} and the whole object is then adjusted accordingly.
#'
#' @examples
#'
#' \dontrun{
#' # Get some data.
#' data("example_mmt")
#'
#' # Let's subset to contig 1, 7675, and 69676.
#' mt1 <- mt_subset(example_mmt,sub_genes = "contig %in% c(1,7675,69676)")
#' mt1
#'
#' # Let's subset to specific organism.
#' mt2 <- mt_subset(example_mmt,sub_samples = "Organism == 'Brocadia'")
#' mt2
#'
#' # Let's do both and remove genes with less than 10000 reads in total.
#' mt3 <- mt_subset(example_mmt,
#'   sub_samples = "Organism == 'Brocadia'",
#'   sub_genes   = "contig %in% c(1,7675,69676)",
#'   minreads    = 10000)
#' mt3
#'
#' # You can also normalise the data and subset.
#' mt4 <- mt_subset(example_mmt,
#'   sub_samples = "Organism == 'Brocadia'",
#'   sub_genes   = "contig %in% c(1,7675,69676)",
#'   minreads    = 10000,
#'   normalise   = "libsize")
#' mt4 # Note "Normalised:" is now included.
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

mt_subset <- function(mmt,sub_genes = NULL,sub_samples = NULL,minreads = 0,normalise = "none"){

  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(mmt$mtmeta) %>% as.numeric()
  ngenesbefore   <- nrow(mmt$mtdata) %>% as.numeric()

  # remove genes below minreads.
  wh <- rowSums(mmt$mtdata[,-1]) >= minreads
  mmt$mtdata <- mmt$mtdata[wh,,drop = F]
  mmt$mtgene  <- if(!is.null(mmt$mtgene)) mmt$mtgene <- mmt$mtgene[wh,,drop = F]

  if(normalise != "none"){
    attributes(mmt)$normalised <- normalise

    # Normalise data.
    if (normalise == "total"){
      mmt$mtdata[,-1] <- as.data.frame(apply(mmt$mtdata[,-1, drop = FALSE],
                          2, function(x) x/sum(x) * 100))
    } else if (normalise == "libsize") {
      mmt$mtdata[,-1] <- mmt$mtdata[,-1] %>%
        as.matrix() %>%
        {t(t(.)/DESeq2::estimateSizeFactorsForMatrix(.))} %>%
        as.data.frame()
    } else {
      stop("normalise: please specify a valid argument.",call. = FALSE)
    }
  }

  # Subset samples.
  if (!is.null(sub_samples)){
    wh <- tryCatch(expr = {
      mmt$mtmeta %>% filter(eval(parse(text = sub_samples)))
    },error = function(e){
      stop("The provided 'sub_samples' string is not meaningfull for the metadata.")
    }) %>% {mmt$mtmeta$SampleID %in% .$SampleID}
    mmt$mtmeta <- mmt$mtmeta[wh,,drop = F]
    mmt$mtdata <- mmt$mtdata[,c(TRUE,wh),drop = F]
  }

  # Subset genes.
  if (!is.null(mmt$mtgene)){
    if(!is.null(sub_genes)){
      wh <- tryCatch(expr = {
        mmt$mtgene %>% filter(eval(parse(text = sub_genes)))
      },error = function(e){
        stop("The provided 'sub_genes' string is not meaningfull for the gene data.")
      }) %>% {mmt$mtgene$GeneID %in% .$GeneID}
      mmt$mtdata <- mmt$mtdata[wh,,drop = F]
      mmt$mtgene <- mmt$mtgene[wh,,drop = F]
    }
  } else {
    stop("There is no gene data available for this mmt object.")
  }

  # Report the removed samples/genes.
  nsamplesafter <- nrow(mmt$mtmeta) %>% as.numeric()
  ngenesafter   <- nrow(mmt$mtdata) %>% as.numeric()

  message(paste(nsamplesbefore - nsamplesafter, "samples and",
                ngenesbefore - ngenesafter, "genes have been filtered \nBefore:",
                nsamplesbefore, "samples and", ngenesbefore, "genes\nAfter:",
                nsamplesafter, "samples and", ngenesafter, "genes"))
  return(mmt)
}
