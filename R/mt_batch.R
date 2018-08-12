#' Correct for batch-effects in \code{mmt} objects based on variable/sample metadata.
#'
#' @title Correct for batch-effects in \code{mmt} objects based on variable/sample metadata.
#'
#' @usage mt_batch(data, ...)
#'
#' @param mmt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param batch (\emph{required}) Variable in metadata to perform batch-correct on.
#' @param method (\emph{optional}) Specify the method for batch-correction. One of:
#' \itemize{
#'    \item \code{"ComBat"}: (\emph{default}) See \link[sva]{ComBat} for details and additional parameter settings.
#' }
#' @param ... Additional parameters for the specified method.
#'
#' @return A modifed mt object
#'
#' @importFrom sva ComBat
#'
#' @export
#'
#' @details The function \code{\link{mt_batch}} tries to eliminate batch-effects.
#'
#' @examples
#'
#' \dontrun{
#' # EXAMPLE GOES HERE
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
mt_batch <- function(mmt,batch,method = "ComBat",...){

  if (!(batch %in% colnames(mmt$mtmeta))) stop("'batch' is not a column in metadata.")
  if (method == "ComBat"){
    mmt$mtdata[,-1] <- sva::ComBat(
      dat   = as.matrix(mmt$mtdata[,-1]),
      batch = mmt$mtmeta[[batch]],...)
  } else {
    stop("Please specify a valid 'method'. See documentation.")
  }

  attributes(mmt)$batch <- list(Var = batch,method = method)

  return(mmt)
}
