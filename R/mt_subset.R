#' Subset mt objects based on sample metadata
#'
#' @title Normalise and subset the data in \code{mmt} objects based on variable/sample metadata.
#'
#' @usage mt_subset(data, ...)
#'
#' @param mmt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param sub_genes (\emph{optional}) A string specifying a logical row subset operation on the mtgene dataframe in the mt object parsed to \link[base]{subset}.
#' @param sub_samples (\emph{optional}) A string specifying a logical row subset operation on the mtmeta dataframe in the mt object parsed to \link[base]{subset}.
#' @param minreads Minimum average number of reads pr. gene. Genes below this value will be removed. (\emph{default:} \code{0})
#' @param frac0 Fraction of zeros allowed per gene. Genes with a higher fraction of zeros will be removed. (\emph{default:} \code{1})
#' @param normalise Normalise the read counts AFTER reads have been removed by the minreads argument but BEFORE any sample/gene subsetting. (\emph{default:} \code{"none"})
#' \itemize{
#'    \item \code{"quantile"}: Quantile normalization. See \link[preprocessCore]{normalize.quantiles} for details.
#'    \item \code{"TPM"}: Normalise read counts to Transcripts Per Milion (TPM).
#'    \item \code{"abundance"}: Normalise read counts to relative abundance.
#'    \item \code{"libsize"}: Normalise the read counts to adjust for gene dispersion and total read counts per sample. See \link[DESeq2]{estimateSizeFactorsForMatrix} for details.
#'    \item \code{"vst"}: Normalise as "libsize" and perform robust log2-transformation. See \link[DESeq2]{vst} for details.
#'    \item \code{"log2"}: Perform log2(x + 1)-transformation.
#'    \item \code{"none"}: No normalisation.
#'    }
#'
#' @return A modifed mt object
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @export
#'
#' @details The function \code{\link{mt_subset}} operates in three steps:
#' \itemize{
#'    \item \strong{Filter} - genes according to \code{minreads} and \code{frac0}.
#'    \item \strong{Normalise} - as specified in \code{normalise}.
#'    \item \strong{Subset} - according to \code{sub_genes} and \code{sub_samples}.
#' }
#'
#' Subsetting are performed on the mtmeta/mtgene data by \link[base]{subset} and the whole object is then adjusted accordingly.
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
mt_subset <- function(mmt,sub_genes = NULL,sub_samples = NULL,minreads = 0,frac0 = 1,normalise = "none"){

  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(mmt$mtmeta) %>% as.numeric()
  ngenesbefore   <- nrow(mmt$mtdata) %>% as.numeric()

  ##### FILTERING #####
  # remove genes below minreads and more than frac0 zeros.
  if (frac0 < 0 | frac0 > 1) stop("'frac0' has to be between [0,1]",call. = FALSE)

  wh <- mmt$mtdata[,.(
    Sum     = Reduce(`+`,.SD), # Sum rows.
    Zeros   = Reduce(`+`,lapply(.SD,`==`,e2 = 0))/length(.SD)), # fraction of zeros per row.
    .SDcols = colnames(mmt$mtdata)[-1]][,Sum >= minreads & Zeros <= frac0]

  mmt$mtdata <- mmt$mtdata[wh]

  if(!is.null(mmt$mtgene)) mmt$mtgene <- mmt$mtgene[wh]

  ##### NORMALISE #####
  if(normalise != "none"){
    if(!is.null(attributes(mmt)$normalised)){
      stop("The data has allready been normalised.",call. = FALSE)
    }
    attributes(mmt)$normalised <- normalise
    Cols <- colnames(mmt$mtdata)[-1]

    # Normalise data.
    if (normalise == "quantile"){
      mmt$mtdata[,(Cols) := data.table(preprocessCore::normalize.quantiles(as.matrix(.SD))),.SDcols = Cols]
    } else if (normalise == "abundance"){
      mmt$mtdata[,(Cols) := lapply(.SD,function(x){ x/sum(x) * 100 }),.SDcols = Cols]
    } else if (normalise == "TPM"){
      if(!("length" %in% colnames(mmt$mtgene))) stop("To normalise by TPM you need a column named 'length' in mtgene, specifying the gene length.")
      TPM <- function(counts,lengths){
        rate = log(counts) - log(as.numeric(lengths))
        exp(rate - log(sum(exp(rate))) + log(10 ^ 6))
      }
      mmt$mtdata[,length := mmt$mtgene$length][,(Cols) := lapply(.SD,TPM,lengths = length),.SDcols = Cols][,length := NULL]
    } else if (normalise == "libsize") {
      mmt$mtdata[,(Cols) := data.table(as.matrix(.SD) %>% {t(t(.)/DESeq2::estimateSizeFactorsForMatrix(.))}),.SDcols = Cols]
    } else if (normalise == "vst"){
      mmt$mtdata[,(Cols) := data.table(as.matrix(.SD) %>% DESeq2::vst()),.SDcols = Cols]
    } else if (normalise == "log2"){
      mmt$mtdata[,(Cols) := lapply(.SD,function(x){log2(x + 1)}),.SDcols = Cols]
    } else {
      stop("normalise: please specify a valid argument.",call. = FALSE)
    }
  }

  ##### SUBSET #####
  # Subset samples.
  if (!is.null(sub_samples)){
    wh <- tryCatch(expr = {
      mmt$mtmeta %>% filter(eval(parse(text = sub_samples)))
    },error = function(e){
      stop("The provided 'sub_samples' string is not meaningfull for the metadata.")
    }) %>% {mmt$mtmeta$SampleID %in% .$SampleID}
    mmt$mtmeta <- mmt$mtmeta[wh]
    mmt$mtdata[,(1+which(!wh)) := NULL]
  }

  # Subset genes.
  if (!is.null(mmt$mtgene)){
    if(!is.null(sub_genes)){
      wh <- tryCatch(expr = {
        mmt$mtgene[eval(parse(text = sub_genes))]
      },error = function(e){
        stop("The provided 'sub_genes' string is not meaningfull for the gene data.")
      }) %>% .$GeneID
      mmt$mtdata <- mmt$mtdata[GeneID %in% wh]
      mmt$mtgene <- mmt$mtgene[GeneID %in% wh]
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
