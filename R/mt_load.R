#' mt_load
#'
#' @title Load the transcriptomics data, associated gene data and metadata.
#'
#' @usage mt_load(mtdata,mtgene = NULL,mtmeta = NULL)
#'
#' @param mtdata (\emph{required}) data.table with read expressions.
#' @param mtgene  (\emph{optional}) data.table with metadata associated to genes (rows).
#' @param mtmeta  (\emph{optional}) data.table with metadata associated to samples (columns).
#'
#' @return A list of class \code{"mt"} with 3 elements.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table setcolorder data.table setnames
#' @importFrom dplyr rename
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' @details If objects of class \code{data.frame} are provided they will be coerced
#' into \code{data.table}. The data structure remains unchanged in most cases,
#' however it is recomended to have the data in \code{data.table} format before
#' loading it.
#'
#' @section The mtdata:
#' The mtdata-table contains read expressions for all genes in each sample. The provided mtdata-table must be a data frame with the following requirements:
#'
#' \itemize{
#'   \item The rows are gene IDs and the columns are samples.
#'   \item The gene ID's are expected to be in eiher the rownames of the data frame or in a column called "ID". Otherwise the function will stop with a message.
#'   \item The column names of the data frame are the sample IDs, exactly matching those in the metadata.
#'   \item Generally avoid special characters and spaces in row- and column names.
#' }
#'
#' A minimal example is available with \code{data("example_mtdata")}.
#'
#' @section The mtgene:
#' The mtgene-table contains additional information about the genes, for example reference database ID's, product, and function, which can be used during analysis. The amount of information in the mtgene-table is unlimited, it can contain any number of columns (variables), however there are a few requirements:
#'
#' \itemize{
#'   \item The gene IDs must be in the first column. These must match exactly to those in the mtdata-table.
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
#'   \item The sample IDs must be in the first column. These sample IDs must match exactly to those in the mtdata-table.
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
#' @examples
#'
#' \dontrun{
#' # Load the different components.
#' data("example_mtmeta")
#' data("example_mtgene")
#' data("example_mtdata")
#'
#' # Combine in one object of class 'mmt'.
#' mt <- mt_load(mtdata = example_mtdata,mtgene = example_mtgene,mtmeta = example_mtmeta)
#'
#' #Show a short summary about the data by simply typing the name of the object in the console
#' mt
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}

mt_load <- function(mtdata,mtgene = NULL,mtmeta = NULL){

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

