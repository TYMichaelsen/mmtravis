#' mt_load
#'
#' This function reads gene expression data, associated annotation and metadata, and returns a list for use in all mmtranscriptome functions. It is therefore required to load data with \code{\link{mt_load}} before any other mmtranscriptome functions can be used.
#'
#' @param mtcount (\emph{required}) data.frame with read counts.
#' @param mtgene  (\emph{optional}) data.frame with metadata associated to genes (rows).
#' @param mtmeta  (\emph{optional}) data.frame with metadata associated to samples (columns).
#' @param mtstat  (\emph{optional}) data.frame with sequencing stats associated to samples (columns).
#'
#' @return A list of class \code{"mt"} with 3 to 4 elements.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' @details The \code{\link{mt_load}} function validates and corrects the provided data frames in different ways to make it suitable for the rest of the mmtranscriptomics workflow. It is important that the provided data frames match the requirements as described in the following sections to work properly.
#'
#' @section The mtcount:
#' The mtcount-table contains read counts for all genes in each sample. The provided mtcount-table must be a data frame with the following requirements:
#'
#' \itemize{
#'   \item The rows are gene IDs and the columns are samples.
#'   \item The gene ID's are expected to be in eiher the rownames of the data frame or in a column called "ID". Otherwise the function will stop with a message.
#'   \item The column names of the data frame are the sample IDs, exactly matching those in the metadata.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' A minimal example is available with \code{data("example_mtcount")}.
#'
#' @section The mtgene:
#' The mtgene-table contains additional information about the genes, for example reference database ID's, product, and function, which can be used during analysis. The amount of information in the mtgene-table is unlimited, it can contain any number of columns (variables), however there are a few requirements:
#'
#' \itemize{
#'   \item The gene IDs must be in the first column. These must match exactly to those in the mtcount-table.
#'   \item Column classes matter, categorical variables should be loaded either \code{as.character()} or \code{as.factor()}, and continuous variables \code{as.numeric()}. See below.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' The \code{\link{mt_load}} function will automatically use the sample IDs in the first column as rownames, but it is important to also have an actual column with sample IDs, so it is possible to fx group by that column during analysis. Any unmatched samples between the otutable and metadata will be removed.
#'
#' A minimal example is available with \code{data("example_mtgene")}.
#'
#' @section The mtmeta:
#' The mtmeta-table contains additional information about the samples, for example where each sample was taken, date, pH, treatment etc, which is used to compare and group the samples during analysis. The amount of information in the mtmeta-table is unlimited, it can contain any number of columns (variables), however there are a few requirements:
#'
#' \itemize{
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the mtcount-table.
#'   \item Column classes matter, categorical variables should be loaded either \code{as.character()} or \code{as.factor()}, and continuous variables \code{as.numeric()}. See below.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' If for example a column is named "Year" and the entries are simply entered as numbers (2011, 2012, 2013 etc), then R will automatically consider these as numerical values (\code{as.numeric()}) and therefore the column as a continuous variable, while it is a categorical variable and should be loaded \code{as.factor()} or \code{as.character()} instead. This has consequences for the analysis as R treats them differently. Therefore either use the \code{colClasses = } argument when loading a csv file or \code{col_types = } when loading an excel file, or manually adjust the column classes afterwards with fx \code{metadata$Year <- as.character(metadata$Year)}.
#'
#' The \code{\link{mt_load}} function will automatically use the sample IDs in the first column as rownames, but it is important to also have an actual column with sample IDs, so it is possible to fx group by that column during analysis. Any unmatched samples between the otutable and metadata will be removed.
#'
#' A minimal example is available with \code{data("example_mtmeta")}.
#'
#' @section The mtstat:
#' The mtstat-table contains sequencing statistics for the samples, and must comply to the following requirements:
#'
#' \itemize{
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the mtcount-table.
#'   \item A column named "Raw", with the raw number of sequences before QC and filtering.
#'   \item A column named "QC", with the number of sequences after QC.
#'   \item A column named "filtered", with the number of sequences after additional removal of rRNA.
#' }
#'
#' The mtstat-table is used for QC assesment.
#'
#' A minimal example is available with \code{data("example_mtstat")}.
#'
#' @examples This
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}

mt_load <- function(mtcount,mtgene = NULL,mtmeta = NULL,mtstat = NULL){

  out <- list(count = NULL,meta = NULL,gene = NULL,stat = NULL)

  ### CHECK INPUT DATA ###
  # Generate dummy meta data if needed.
  if (is.null(mtmeta)){
    mtmeta <- data.frame(
      SampleID = colnames(mtcount[,-1]),
      Dummy    = "All")
    warning("No sample metadata provided, creating dummy metadata.\n", call. = FALSE)
  }

  # Check gene data.
  gflag <- TRUE
  if (is.null(mtgene)){
    gflag  <- FALSE
    warning("No gene data provided.\n", call. = FALSE)
  }

  # Check stat data.
  sflag <- TRUE
  if (is.null(mtstat)){
    sflag  <- FALSE
    warning("No stat data provided.\n", call. = FALSE)
  }

  ### CORRECT NAMING ###
  colnames(mtcount) <- stringr::str_replace_all(colnames(mtcount), "[^[:alnum:]]", "_")
  colnames(mtmeta)  <- stringr::str_replace_all(colnames(mtmeta), "[^[:alnum:]]", "_")
  mtmeta[,1]        <- stringr::str_replace_all(mtmeta[,1], "[^[:alnum:]]", "_")
  if(gflag) colnames(mtgene) <- stringr::str_replace_all(colnames(mtgene), "[^[:alnum:]]", "_")
  if(sflag){
    colnames(mtstat) <- stringr::str_replace_all(colnames(mtstat), "[^[:alnum:]]", "_")
    mtstat[,1]       <- stringr::str_replace_all(mtstat[,1], "[^[:alnum:]]", "_")
  }

  ### ENSURE CORRECT NAMING ESSENTIAL COLUMNS ###
  # mtcount
  i <- which(colnames(mtcount) == "GeneID")[1]
  if (!is.na(i) & i != 1){
    mtmeta <- dplyr::rename(mtcount,GeneID_1 = GeneID)
    colnames(mtcount)[1] <- "GeneID"
    message("You had a column named 'GeneID' in mtcount which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
  } else {
    colnames(mtcount)[1] <- "GeneID"
  }
  # mtmeta
  i <- which(colnames(mtmeta) == "SampleID")[1]
  if (!is.na(i) & i != 1){
    mtmeta <- dplyr::rename(mtmeta,SampleID_1 = SampleID)
    colnames(mtmeta)[1] <- "SampleID"
    message("You had a column named 'SampleID' in mtmeta which wasn't the first column. Renaming it to 'SampleID_1' to avoid conflicts.")
  } else {
    colnames(mtmeta)[1] <- "SampleID"
  }
  # mtstat
  if (sflag){
    if (any(!c("Raw","QC","filtered") %in% colnames(mtstat))) stop("'mtstat' are missing vital columns or name is wrong. Please read documentation.")
    i <- which(colnames(mtstat) == "SampleID")[1]
    if (!is.na(i) & i != 1){
      mtstat <- dplyr::rename(mtstat,SampleID_1 = SampleID)
      colnames(mtstat)[1] <- "SampleID"
      message("You had a column named 'SampleID' in mtstat which wasn't the first column. Renaming it to 'SampleID_1' to avoid conflicts.")
    } else {
      colnames(mtstat)[1] <- "SampleID"
    }
  }
  # mtgene
  if (gflag){
    i <- which(colnames(mtgene) == "GeneID")[1]
    if (!is.na(i) & i != 1){
      mtstat <- dplyr::rename(mtstat,GeneID_1 = GeneID)
      colnames(mtgene)[1] <- "GeneID"
      message("You had a column named 'GeneID' in mtgene which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
    } else {
      colnames(mtgene)[1] <- "GeneID"
    }
  }

  ### CHECK THAT THINGS WENT WELL ###
  # Check consistent sample names.
  i <- list(colnames(mtcount[,-1]),mtmeta[,1],if(sflag) mtstat[,1] else NULL) %>%
    {.[!sapply(.,is.null)]} %>% outer(.,.,Vectorize(setequal)) %>% all()
  if ( !i ){
    stop("Sample names are not matching 1:1 or misspecified in 'mtcount', 'mtstat', and/or 'mtmeta'. Please read the documentation.")
  }
  # Check consistent gene names.
  i <- if(gflag) setequal(mtgene[,1],mtcount[,1]) else TRUE
  if ( !i ){
    stop("Gene names are not matching 1:1 or misspecified in 'mtcount' and/or 'mtgene'. Please read the documentation.")
  }

  ### DUMP OUTPUT ###
  sampOrd <- as.character(mtmeta$SampleID)

  out$count <- mtcount[,c(1,1+match(sampOrd,colnames(mtcount)[-1]))]
  out$meta  <- mtmeta[match(sampOrd,mtmeta$SampleID),]
  out$stat  <- if(sflag) mtstat[match(sampOrd,mtstat$SampleID),] else NULL
  out$gene  <- if(gflag) mtgene[match(mtcount$GeneID,mtgene$GeneID),] else NULL


  return(out)
}
