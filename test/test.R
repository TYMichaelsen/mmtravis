# TEST SCRIPT
library(magrittr)
library(readxl)
library(tidyverse)
library(mmtravis)
library(data.table)


tab <- fread(
  file             = "test/counts.txt",
  header           = T,
  stringsAsFactors = F,
  check.names      = F)
metaVars <- c("Genome","Contig","start","end","strand","contig_len","ftype","gene","EC_number","product","locus_tag","function","inference")

# Prepare the gene data.
genedata <- select(tab,1,2) %>%
  separate(.,
           col    = 2,
           into   = metaVars,
           sep    = "[|]",
           remove = T) %>%
  as.data.frame() %>%
  mutate(length = as.numeric(end) - as.numeric(start))

# Prepare the count table.
counttab <- select(tab,-2) %>%
  as.data.frame()

# Metadata.
metadata  <- read_excel("test/metadata.xlsx") %>%
  as.data.frame() %>%
  select(SeqID,everything())

### Load data.
mt <- mt_loadMetaT(counts.txt = "test/counts.txt",seqstat.txt = "test/seqstat.txt",mtmeta = mtmeta)

mt <- mt_load2(mtdata = counttab,mtgene = genedata,mtmeta = metadata)
### Test functions. ############################################################
# Subset.
de <- mt_subset2(mt,minreads = 0,frac0 = 1,normalise = "libsize")

mt_subset(mt,sub_genes = "Contig == '1'")

mt_plotpairs(mt,samples = c("HQ180523_13","HQ180523_14"),label_by = "Label",linesize = 2)

hest <- mt_gather(de)

de_bch <- mt_batch(de,batch = "VFA")


mt_stackbar(mt,detailed = T)

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
  group      = "Type",
  row_labels = c("GeneID","product"),
  intercept  = "ANAMMOX")

# Load data.
data("example_mmt")

# Compute the statistics.
DE_mt <- mt_diffexprs(example_mmt,
                      group      = "Type",
                      row_labels = c("GeneID","product"),
                      intercept  = "ANAMMOX")

# Extract results for given levels.
mt_res <- mt_dumpDE(DE_mt,nom = "ELECTRODE",denom = "SUSPENTION")

# Show the output.
mt_res$MAplot
mt_res$BOXplot
head(mt_res$Table)




