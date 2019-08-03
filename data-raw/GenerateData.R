library(data.table)
library(magrittr)

################################################################################
### Code to generate the COG database ##########################################
################################################################################

# FILES USED FROM FOLDER /COG_data/:
# cognames2003-2014-tab   Downloaded from NCBI FTP site ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data on d. 02-08-2019
# fun2003-2014.tab        Downloaded from NCBI FTP site ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data on d. 02-08-2019

create <- function(map,description){
  out <- merge(map,description,by = colnames(map)[1],all.x = T) %>%
    plyr::dlply(.variables = colnames(map)[1],.fun = function(dt){list(COGs = dt$COG,Description = unique(dt$Description))})
  attr(out,"split_type")   <- NULL
  attr(out,"split_labels") <- NULL
  keep <- lapply(out,"[[","COGs") %>% {sapply(.,length) > 0}
  return(out[keep])
}

cog_to_code <- fread("data-raw/COG_data/cognames2003-2014.tab",col.names = c("COG","Code1","COG_desc")) %>%
  .[,c("X1","X2","X3","X4") := tstrsplit(Code1,"",fill = NA)] %>%
  melt(id.vars = 1:3,measure.vars = 4:7,na.rm = T,value.name = "Code") %>%
  .[,c("Code","COG")]

code_nm <- fread("data-raw/COG_data/fun2003-2014.tab",col.names = c("Code","Description"))

COG_DB <- list(Category = create(cog_to_code,code_nm),
               AllCOG   = unique(cog_to_code$COG))

# Remove category S and R, they are quite meaningless.
COG_DB$Category$S <- NULL
COG_DB$Category$R <- NULL

# Add category for unannotated.
COG_DB$Category$None <- list(COGs = "None",Description = "Not annotated")

################################################################################
### Code to generate the KO database ###########################################
################################################################################

# FILES USED FROM FOLDER /KO_data/:
# *_to_ko.tsv   Files mapping from identifiers to KOs. Created using the script "commands.txt".
# *_names.tsv   Files mapping from identifiers to descriptions. Created using the script "commands.txt".
# all_ko.tsv    File containing a list of all KOs. Created using the script "commands.txt".

create <- function(map,description){
  out <- merge(map,description,by = colnames(map)[1],all.x = T) %>%
    plyr::dlply(.variables = colnames(map)[1],.fun = function(dt){list(KOs = dt$KO,Description = unique(dt$Description))})
  attr(out,"split_type")   <- NULL
  attr(out,"split_labels") <- NULL
  keep <- lapply(out,"[[","KOs") %>% {sapply(.,length) > 0}
  return(out[keep])
}

# Load the data.
all_KO         <- fread("data-raw/KO_data/all_ko.tsv",sep = "\t")

path_to_ko     <- fread("data-raw/KO_data/pathway_to_ko.tsv",       sep = "\t")
module_to_ko   <- fread("data-raw/KO_data/module_to_ko.tsv",        sep = "\t")
rclass_to_ko   <- fread("data-raw/KO_data/reaction-class_to_ko.tsv",sep = "\t")
reaction_to_ko <- fread("data-raw/KO_data/reaction_to_ko.tsv",      sep = "\t")
EC_to_ko       <- fread("data-raw/KO_data/ec_to_ko.tsv",            sep = "\t")

path_nm     <- fread("data-raw/KO_data/pathway_names.tsv",       sep = "\t")
module_nm   <- fread("data-raw/KO_data/module_names.tsv",        sep = "\t")
rclass_nm   <- fread("data-raw/KO_data/reaction-class_names.tsv",sep = "\t")
reaction_nm <- fread("data-raw/KO_data/reaction_names.tsv",      sep = "\t")
ec_nm       <- fread("data-raw/KO_data/ec_names.tsv",            sep = "\t")

KEGG_DB <- list(
  pathway        = create(path_to_ko,path_nm),
  module         = create(module_to_ko,module_nm),
  reaction_class = create(rclass_to_ko,rclass_nm),
  reaction       = create(reaction_to_ko,reaction_nm))

# Add category for unannotated.
KEGG_DB <- lapply(KEGG_DB,`c`,list(None = list(KOs = "None",Description = "Not annotated")))

# Add list with all KOs.
KEGG_DB$AllKO = all_KO$KO

### SAVE ALL DATA ###

usethis::use_data(KEGG_DB,COG_DB,internal = T,overwrite = T)
