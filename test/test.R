# TEST SCRIPT
library(magrittr)
library(readxl)
library(tidyverse)
library(mmtranscriptome)

### Load data.
mtmeta  <- read_excel("test/metadata_2.xlsx") %>%
  as.data.frame() %>%
  mutate(Samplename = gsub("_SP","",Samplename)) %>%
  separate(
    col  = Samplename,
    into = c("Genus","Replicate","Type"),
    sep  = "_")

mtstat <- read.delim(
  file             = "test/seqstat_2.txt",
  header           = T,
  stringsAsFactors = F)

count_tab <- read.delim(
  file             = "test/counts_2.txt",
  header           = T,
  stringsAsFactors = F,
  check.names      = F)
metaVars <- colnames(count_tab)[2] %>%
  gsub("[{]|[}]","",x = .) %>%
  {strsplit(.,split = ";")[[1]]}

mtgene <- select(count_tab,1,2) %>%
  separate(.,
    col    = 2,
    into   = metaVars,
    sep    = "[|];[|]",
    remove = T)

mtcount <- select(count_tab,-2)

### Test functions.

mt <- mt_load(mtcount = mtcount,mtgene = mtgene,mtstat = mtstat,mtmeta = mtmeta)


mt_stats(mt,group_by = "Genus")


mt_binstats <- function(mt){

  out <- list()

  ### GET SOME BIN STATS, IF POSSIBLE ###
  if(!is.null(mt$gene) & any(tolower(colnames(mt$gene)) %in% c("contig","scaffold"))){
    wh <- grep("contig|scaffold",colnames(mt$gene),value = T,ignore.case=TRUE)[1]

    # Get number of contigs/scaffolds.
    out$no_bins <- mt$gene %>%
      .[[wh]] %>%
      unique() %>%
      length()
    # Get average gene expression of each bin.
    avg <- gather(mt$count,key = SampleID,value = Expr,-GeneID) %>%
      left_join(.,mt$gene,by = "GeneID") %>%
      group_by_(wh,"SampleID") %>%
      summarize(avg = median(Expr)) %>%
      filter(avg != 0)

    out$expr_bins

    de <- avg %>% group_by(SampleID) %>% summarize(avg_sample = mean(avg))
    ggplot(de,aes(x = SampleID,y = avg)) + geom_boxplot()

  }

  ### GET SOME SEQUENCING STATS, IF POSSIBLE ###
  if(!is.null())

}


mt_stats(mt,group_by = "Genus")

