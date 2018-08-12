
mt_subset2 <- function(mmt,sub_genes = NULL,sub_samples = NULL,minreads = 0,frac0 = 1,normalise = "none"){

  #For printing removed samples and OTUs
  nsamplesbefore <- nrow(mmt$mtmeta) %>% as.numeric()
  ngenesbefore   <- nrow(mmt$mtdata) %>% as.numeric()

  # Copy to utilize by reference.
  mmt <- lapply(mt,copy) %>% `class<-`("mmt")

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
    } else if (normalise == "total"){
      mmt$mtdata[,(Cols) := lapply(.SD,function(x){ x/sum(x) * 100 }),.SDcols = Cols]
    } else if (normalise == "TPM"){
      if(!("length" %in% colnames(mmt$mtgene))) stop("To normalise by TPM you need a column named 'length' in mtgene, specifying the gene length.")
      TPM <- function(counts,lengths){
        rate = log(counts) - log(lengths)
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
