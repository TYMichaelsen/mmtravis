
mt_subset2 <- function(mmt,sub_genes = NULL,sub_samples = NULL,minreads = 0,frac0 = 1,normalise = "none"){

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
    }) %>% .$SampleID
    mmt$mtmeta <- mmt$mtmeta[SampleID %in% wh]
    mmt$mtdata[,(wh) := NULL]
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

################################################################################
################################################################################
################################################################################
mt_load2 <- function(mtdata,mtgene = NULL,mtmeta = NULL){

  ### CHECK INPUT DATA ###
  # Generate meta data if needed.
  if (is.null(mtmeta)){
    mtmeta <- data.table(
      SampleID = colnames(mtdata[,-1]))
    warning("No sample metadata provided, creating dummy metadata.\n", call. = FALSE)
  } else {
    mtmeta <- data.table(mtmeta)
  }

  # Generate gene data if needed.
  if (is.null(mtgene)){
    mtgene <- data.table(
      GeneID = mtdata[,1])
    warning("No gene data provided, creating dummy genedata.\n", call. = FALSE)
  } else {
    mtgene <- data.table(mtgene)
  }

  # Set mtdata to correct class.
  mtdata <- data.table(mtdata)

  ### CORRECT NAMING ###
  setnames(mtdata,old = colnames(mtdata),new = stringr::str_replace_all(colnames(mtdata), "[^[:alnum:]]", "_"))
  setnames(mtmeta,old = colnames(mtmeta),new = stringr::str_replace_all(colnames(mtmeta), "[^[:alnum:]]", "_"))
  setnames(mtgene,old = colnames(mtgene),new = stringr::str_replace_all(colnames(mtgene), "[^[:alnum:]]", "_"))
  mtmeta[,1] <- stringr::str_replace_all(mtmeta[[1]], "[^[:alnum:]]", "_")

  ### ENSURE CORRECT NAMING ESSENTIAL COLUMNS ###
  # mtdata
  i <- which(colnames(mtdata) == "GeneID")[1]
  if (!is.na(i) & i != 1){
    setnames(mtdata,old = "GeneID",new = "GeneID_1")
    setnames(mtdata,old = 1,new = "GeneID")
    message("You had a column named 'GeneID' in mtdata which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
  } else {
    setnames(mtdata,old = 1,new = "GeneID")
  }
  # mtmeta
  i <- which(colnames(mtmeta) == "SampleID")[1]
  if (!is.na(i) & i != 1){
    setnames(mtmeta,old = "SampleID",new = "SampleID_1")
    setnames(mtmeta,old = 1,new = "SampleID")
    message("You had a column named 'SampleID' in mtmeta which wasn't the first column. Renaming it to 'SampleID_1' to avoid conflicts.")
  } else {
    setnames(mtmeta,old = 1,new = "SampleID")
  }
  # mtgene
  i <- which(colnames(mtgene) == "GeneID")[1]
  if (!is.na(i) & i != 1){
    setnames(mtgene,old = "GeneID",new = "GeneID_1")
    setnames(mtgene,old = 1,new = "GeneID")
    message("You had a column named 'GeneID' in mtgene which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
  } else {
    setnames(mtgene,old = 1,new = "GeneID")
  }

  ### CHECK THAT THINGS WENT WELL ###
  # Check consistent sample names.
  i <- setequal(colnames(mtdata[,-1,drop = F]),mtmeta[[1]])
  if ( !i ){
    stop("Sample names are not matching 1:1 or misspecified in 'mtdata' and/or 'mtmeta'. Please read the documentation.")
  }
  # Check consistent gene names.
  i <- setequal(mtgene[[1]],mtdata[[1]])
  if ( !i ){
    stop("Gene names are not matching 1:1 or misspecified in 'mtdata' and/or 'mtgene'. Please read the documentation.")
  }

  ### DUMP OUTPUT ###
  sampOrd <- as.character(mtmeta$SampleID)

  setcolorder(mtdata,neworder = c(1,1+match(sampOrd,colnames(mtdata)[-1])))
  setroworder(mtmeta,neworder = match(sampOrd,mtmeta$SampleID))
  setroworder(mtgene,neworder = match(mtdata$GeneID,mtgene$GeneID))

  out <- list(
    mtdata = mtdata,
    mtgene = mtgene,
    mtmeta = mtmeta)

  class(out) <- "mmt"
  return(out)
}
