# TEST SCRIPT
library(magrittr)
library(readxl)
library(tidyverse)
library(mmtravis)

### Load data.
mtmeta  <- read_excel("test/metadata.xlsx") %>%
  as.data.frame() %>%
  mutate(Samplename = gsub("_SP","",Samplename)) %>%
  separate(
    col  = Samplename,
    into = c("Genus","Replicate","Type"),
    sep  = "_")

mt <- mt_loadMetaT(counts.txt = "test/counts.txt",seqstat.txt = "test/seqstat.txt",mtmeta = mtmeta)

### Test functions. ############################################################
# Subset.
mt     <- mt_subset(example_mmt,minreads = 5000,normalise = "total")

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




