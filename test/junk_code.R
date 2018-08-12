## Junk code that may be handy later.


# Calculate summary stats.
gnomes <- unique(obj$mtgene$Genome) %w/o% "assembly" %>% length

summar <- data.frame(
  Samples        = ncol(obj$mtdata) - 1,
  Genes          = nrow(obj$mtdata),
  Genomes        = unique(obj$mtgene$Genome) %w/o% "assembly" %>% length,
  'Genes/Genome' =
)


#'    \item \code{"vst"}: First normalise as \code{libsize}, then normalise to variance-stabilized log2 values. See \link[DESeq2]{varianceStabilizingTransformation} for details.
#'    \item \code{"log2"}: First normalise as \code{libsize}, before applying a simple log2(x + 1) transformation.




mt_subset2 <- function(mmt,sub_genes = NULL,sub_samples = NULL,minreads = 0,normalise = "none"){

  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(mmt$mtmeta) %>% as.numeric()
  ngenesbefore   <- nrow(mmt$mtdata) %>% as.numeric()

  # remove genes below minreads.
  wh <- data.table(mmt$mtdata[,-1]) %>%
    .[, `:=`(i = rowMeans(.SD, na.rm=T) > minreads)] %>% .[,i]

  mmt$mtdata <- mmt$mtdata[wh,,drop = F]
  mmt$mtgene  <- if(!is.null(mmt$mtgene)) mmt$mtgene <- mmt$mtgene[wh,,drop = F]

  ##### NORMALIZATION #####
  if(normalise != "none"){
    attributes(mmt)$normalised <- normalise

    # Normalise data.
    if (normalise == "total"){
      mmt$mtdata[,-1] <- as.data.frame(apply(mmt$mtdata[,-1, drop = FALSE],
                                             2, function(x) x/sum(x) * 100))
    } else if (normalise %in% c("libsize","log2")) {
      mmt$mtdata[,-1] <- mmt$mtdata[,-1] %>%
        as.matrix() %>%
        {t(t(.)/DESeq2::estimateSizeFactorsForMatrix(.))} %>%
        as.data.frame()
      if (normalise == "log2"){
        mmt$mtdata[,-1] <- log2(mmt$mtdata[,-1] + 1)
      }
    } else if (normalise == "vst"){
      mmt$mtdata[,-1] <- mmt$mtdata[,-1] %>%
        as.matrix() %>%
        varianceStabilizingTransformation() %>%
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
