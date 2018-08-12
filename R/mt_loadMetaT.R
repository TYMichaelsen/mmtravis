#' mt_loadMetaT
#'
#' @title Load data produced by the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline into R.
#'
#' @usage mt_loadMetaT(counts.txt,seqstat.txt,mtmeta = NULL)
#'
#' @param counts.txt (\emph{required}) The count table as outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline
#' @param seqstat.txt (\emph{required}) The seqstat table as outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline
#' @param mtmeta (\emph{optional}) data.frame with metadata associated to samples (columns). See \code{\link{mt_load}} for details.
#'
#' @return A \code{mmt} object.
#'
#' @details This function is a wrapper to transform .txt files outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline into R.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom dplyr select left_join rename mutate
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
mt_loadMetaT <- function(counts.txt,seqstat.txt,mtmeta = NULL){

# Load the count + genedata.
tab <- fread(
  file             = counts.txt,
  header           = T,
  stringsAsFactors = F,
  check.names      = F)
metaVars <- c("Genome","Contig","start","end","strand","contig_len","ftype","gene","EC_number","product","locus_tag","function","inference")

# Prepare the gene data.
mtgene <- select(tab,1,2) %>%
  separate(.,
           col    = 2,
           into   = metaVars,
           sep    = "[|]",
           remove = T) %>%
  as.data.frame()

# Prepare the count table.
mtdata <- select(tab,-2) %>%
  as.data.frame()

# Load as mmt object.
out <- mt_load(mtdata = mtdata,mtgene = mtgene,mtmeta = mtmeta)

# Add additional sequencing data.
seqstats <- read.delim(
  file             = seqstat.txt,
  header           = T,
  stringsAsFactors = F) %>%
  {`colnames<-`(.,stringr::str_replace_all(colnames(.), "[^[:alnum:]]", "_"))} %>%
  mutate(SeqID = stringr::str_replace_all(SeqID, "[^[:alnum:]]", "_")) %>%
  dplyr::rename(SampleID = SeqID)

out$mtmeta <- left_join(out$mtmeta,seqstats,by = "SampleID")
return(out)
}
