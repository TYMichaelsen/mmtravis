#' mt_stattab
#'
#' This function makes a table of the 'mtstat' component.
#'
#' @param mt       (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param add_meta (\emph{optional}) Add an aditional column of metadata to the table, for explanatory reasons.
#' @param textsize (\emph{optional}) Adjust the textsize of the table.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate_at left_join rename
#' @importFrom tibble rownames_to_column add_column
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.table ttheme_default
#'
#' @export
#'
#' @details The \code{\link{mt_stattab}} function takes a 'mt' object and generates a table of sequencing statistics based on the sequencing stats stored in the object, if such are available. Throws an error if data is missing.
#'
#' @examples Some example.
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}

mt_stattab <- function(mt,add_meta = NULL,textsize = 7){

  if(!is.null(mt$stat)){
    # Create output table.
    mapped <- colSums(mt$count[,-1]) %>%
      data.frame(mapped = .) %>%
      rownames_to_column("SampleID")
    tab <- left_join(mt$stat,mapped,by = "SampleID") %>%
      mutate_at(.vars = c("Raw","QC","filtered","mapped"),.funs = function(x){round(x/10^6,2)}) %>%
    rename("# raw reads" = Raw,
           "# reads \n after QC" = QC,
           "# reads \nrRNA removed" = filtered,
           "# reads \n mapped" = mapped)
  } else {
    stop("You do not have a 'mtstat' object assigned to your data.",call. = FALSE)
  }

  # Group by grouping variable.
  if(!is.null(add_meta)){
    if (all(add_meta %in% colnames(mt$meta))){
      tab <- left_join(mt$meta[,c("SampleID",add_meta)],tab,by = "SampleID")
    } else {
      stop("Your 'add_meta' is not in metadata.",call. = FALSE)
    }
  }

  # Plot the table.
  tt2 <- ttheme_default(
    core      = list(fg_params = list(hjust = 1, x = 0.95, fontsize = textsize)),
    colhead   = list(fg_params = list(fontsize = textsize+1)))

  grid.table(tab, rows = NULL, theme = tt2)
}
