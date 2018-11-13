#' Calculate the phi-statistic for all genes in mmt object.
#'
#' @title Calculates the phi-statistic of all genes in \code{mmt} object.
#'
#' @usage mt_subset(data, ...)
#'
#' @param mmt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#'
#' @return A genes x genes matrix with phi-statistics.
#'
#' @importFrom propr propr
#' @importFrom tibble column_to_rownames
#'
#' @export
#'
#' @details Some details here.
#'
#' @examples
#'
#' \dontrun{
#' # An example goes here.
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
mt_phi <- function(mmt){

  # Check if data is proportional.
  if (!(attributes(mmt)$normalised %in% c("TPM","abundance")))
    stop("Data has to be relative (normalised by 'abundance' or 'TPM')",
         call. = FALSE)

  # Build the matrices holding the data.
  dat <- mmt$mtdata %>%
    as.data.frame() %>%
    column_to_rownames("GeneID") %>%
    as.matrix()

  # Compute phi.
  phi <- propr(counts = t(dat),metric = "phi",symmetrize = TRUE)

  return(phi@matrix)
}
