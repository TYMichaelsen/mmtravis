
library(data.table)

### TRY TO SEARCH KOs IN R ###


### BUILD DATABASE ###

DB <- fread("test/testKO/Module_steps.tsv",
  col.names = c("Module","Steps"),
  header    = F,
  sep       = "\t")

# Make a list of each module, with an element for each step.
DB[,Steps := lapply(Steps,function(x){
  strsplit(x,split = ";")[[1]] %>%
  sapply(strsplit,split = ",") %>%
  {setNames(.,nm = 1:length(.))}
})]

# Add a column with unique KOs for each module.
DB[,Unique := lapply(Steps,function(x){
  unlist(x) %>% unique()
})]


### TESTING ###
KOlist <- readRDS("test/testKO/KOlist_bins.rds")

# Convert to data.frame
KOdf <- sapply(KOlist,paste0,collapse = ";") %>%
  data.frame(KOs = .) %>%
  rownames_to_column("Bin") %>%
  mutate(Bin = as.character(Bin),KOs = as.character(KOs))


quer <- query_genomes_to_modules(
  GENOME_INFO   = KOdf[1,],
  GENOME_ID_COL = "Bin",
  GENES_COL     = "KOs",
  splitBy       = ";")



data(data_example_multi_ECs_KOs) # load example data set

hest <- data_example_multi_ECs_KOs[,c("ID","KOs")]
hest$ID <- paste0("hyrbin",1:5)
names(data_example_multi_ECs_KOs)
# "ID"       "ORG_ID"   "ORGANISM" "KOs"      "ECs"
OUT <- query_genomes_to_modules(hest,GENOME_ID_COL = "ID",
                                GENES_COL = "KOs",MODULE_ID = paste("M0000",1:5,sep=""),
                                META_OUT = T,ADD_OUT = T)


