# TEST SCRIPT
library(magrittr)
library(readxl)
library(tidyverse)
library(mmtravis)
library(data.table)


tab <- fread(
  file             = "test/counts_CP361.txt",
  header           = T,
  stringsAsFactors = F,
  check.names      = F)
metaVars <- c("Genome","Contig","Gene_start","Strand","Inference","Contig_len","ftype","length","Gene","EC_number","COG","product","locus_tag","function")

# Prepare the gene data.
genedata <- select(tab,1,2) %>%
  separate(.,
           col    = 2,
           into   = metaVars,
           sep    = "[|]",
           remove = T) %>%
  as.data.frame()

# Prepare the count table.
counttab <- select(tab,-2) %>%
  as.data.frame()

# Metadata.
metadata  <- read_excel("test/CP361_metadata.xlsx") %>%
  mutate(Samplename = paste0("S",Samplename)) %>%
  mutate(SEQID_2    = ifelse(is.na(SEQID_2),"",SEQID_2))

# As we have two batches, we need to reshape metadata.
bch1 <- metadata %>% select(-SEQID_2) %>% dplyr::rename(SEQID = SEQID_1) %>% mutate(batch = "1")
bch2 <- subset(metadata,SEQID_2 != "") %>% select(-SEQID_1) %>% dplyr::rename(SEQID = SEQID_2) %>% mutate(batch = "2")
metadata2 <- rbind(bch1,bch2)

### Load data.
mt <- mt_loadMetaT(
  counts.txt  = "test/counts_CP361.txt",
  seqstat.txt = "test/seqstat_CP361.txt",
  mtmeta      = metadata2)

mt <- mt_load(mtdata = counttab,mtgene = genedata,mtmeta = metadata2)
### Test functions. ############################################################
# Subset.
wh <- subset(metadata2,batch == 2)$Samplename %>%
  paste0(.,collapse = "','") %>%
  paste0("Samplename %in% c('",.,"')")

batch <- mt_subset(mt,sub_samples = wh)


gat <- mt_gather(mt,metavars = "MetaTrans")

de <- mt_subset2(mt,minreads = 20,frac0 = 0.5,normalise = "TPM",
                 sub_genes = "Contig == '1'",
                 sub_samples = "VFA == 'Acetate'")

de <- mt_subset2(mt,sub_genes = "Contig == '1'")

mt_plotpairs(mt,samples = c("HQ180523_1","HQ180523_2"),label_by = "Label",linesize = 2)

hest <- mt_gather(mt)

de_bch <- mt_batch(de,batch = "VFA")


de <- mt_stackbar(mt,detailed = T,group_by = "MetaTrans")

data("example_mmt")

# Plot all samples (might take some time)
mt_plotpairs(example_mmt)

# Plot replicates to see correlation.
mt_plotpairs(example_mmt,samples = c("HQ180323_13","HQ180323_14"),label_by = "Replicate")

# Ordination.
vis_ordinate(mt,
  sample_color_by = "Genus",
  sample_shape_by = "Replicate")

# Boxplotting.
# Let's plot some random genes.
wh <- example_mmt$mtdata$GeneID[c(1,22,11,53)]


vis_boxplot(example_mmt,
  group_by    = "Type",
  boxgroup_by = "Organism",
  row_show    = 10,
  row_labels  = "product",
  title_newline = T)

vis_boxplot(AalborgWWTPs,group_by = "Period")

# Differential expression stuff.
DE_mt <- mt_diffexprs(mt,
  group      = "Experiment",
  row_labels = c("GeneID","product"),
  intercept  = "mAlg")


mt_res <- mt_dumpDE(DE_mt,num = "mPec",denom = "mAlg",title_size = 10)

mt_res$BOXplot + scale_y_log10()
mt_res$MAplot

# Load data.
data("example_mmt")

# Compute the statistics.
DE_mt <- mt_diffexprs(example_mmt,
                      group      = "Type",
                      row_labels = c("GeneID","product"),
                      intercept  = "ANAMMOX")

# Extract results for given levels.
mt_res <- mt_dumpDE(DE_mt,num = "ELECTRODE",denom = "SUSPENTION",title_size = 8)

# Show the output.
mt_res$MAplot
mt_res$BOXplot
head(mt_res$Table)

# Enrichment analysis.

# Suppose we have a set of genes which want to test for functional enrichment.
Genes <- sample(mt$mtgene$GeneID,size = 1000)
head(Genes)

COGenrich <- mt_enrichCOG(mt,
  GeneIDs     = Genes,
  COGs        = "COG",
  alternative = "two.sided",
  show_p      = 1,
  size_p      = 4,
  unannotated = F)

COGenrich$plot




CP460 <- readRDS("../../../../DNASense ApS Dropbox/DNASense/CP/cp-ongoing/CP460 [BASF Metatranscriptomics]/analysis/output/mt_object.rds")

CP460_CDS <- mt_subset(CP460,sub_genes = "ftype != 'intergenic_region'")

Genes <- mt_subset(CP460,sub_genes = "Genome == 'Bin1'")$mtgene$GeneID

KOenrich <- mt_enrichKO(CP460,
  GeneIDs     = Genes,
  type        = "reaction",
  KOs         = "KO",
  alternative = "two.sided",
  size_p      = 4,
  show_p      = 1,
  unannotated = F)

KOenrich$plot
