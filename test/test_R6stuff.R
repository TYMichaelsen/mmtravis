mmt <- R6::R6Class("mmt",
                   private = list(
                     ..mtcount = NULL,
                     ..mtgene  = NULL,
                     ..mtmeta  = NULL),
                   public = list(
                     ### WRITE THE INITIALIZE FUNCTION ----------------------------------------------
                     initialize = function(mtmeta = NULL,mtcount = NULL,mtgene = NULL){

                       ### CHECK INPUT DATA ###
                       # Generate meta data if needed.
                       if (is.null(mtmeta)){
                         mtmeta <- data.frame(
                           SampleID = colnames(mtcount[,-1]))
                         warning("No sample metadata provided, creating dummy metadata.\n", call. = FALSE)
                       }
                       # Generate gene data if needed.
                       if (is.null(mtgene)){
                         mtgene <- data.frame(
                           GeneID = mtcount[,1])
                         warning("No gene data provided, creating dummy genedata.\n", call. = FALSE)
                       }

                       ### CORRECT NAMING ###
                       colnames(mtcount) <- stringr::str_replace_all(colnames(mtcount), "[^[:alnum:]]", "_")
                       colnames(mtmeta)  <- stringr::str_replace_all(colnames(mtmeta), "[^[:alnum:]]", "_")
                       colnames(mtgene)  <- stringr::str_replace_all(colnames(mtgene), "[^[:alnum:]]", "_")
                       mtmeta[,1]        <- stringr::str_replace_all(mtmeta[,1], "[^[:alnum:]]", "_")

                       ### ENSURE CORRECT NAMING ESSENTIAL COLUMNS ###
                       # mtcount
                       i <- which(colnames(mtcount) == "GeneID")[1]
                       if (!is.na(i) & i != 1){
                         mtmeta <- dplyr::rename(mtcount,GeneID_1 = GeneID)
                         colnames(mtcount)[1] <- "GeneID"
                         message("You had a column named 'GeneID' in mtcount which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
                       } else {
                         colnames(mtcount)[1] <- "GeneID"
                       }
                       # mtmeta
                       i <- which(colnames(mtmeta) == "SampleID")[1]
                       if (!is.na(i) & i != 1){
                         mtmeta <- dplyr::rename(mtmeta,SampleID_1 = SampleID)
                         colnames(mtmeta)[1] <- "SampleID"
                         message("You had a column named 'SampleID' in mtmeta which wasn't the first column. Renaming it to 'SampleID_1' to avoid conflicts.")
                       } else {
                         colnames(mtmeta)[1] <- "SampleID"
                       }
                       # mtgene
                       i <- which(colnames(mtgene) == "GeneID")[1]
                       if (!is.na(i) & i != 1){
                         colnames(mtgene)[1] <- "GeneID"
                         message("You had a column named 'GeneID' in mtgene which wasn't the first column. Renaming it to 'GeneID_1' to avoid conflicts.")
                       } else {
                         colnames(mtgene)[1] <- "GeneID"
                       }

                       ### CHECK THAT THINGS WENT WELL ###
                       # Check consistent sample names.
                       i <- setequal(colnames(mtcount[,-1]),mtmeta[,1])
                       if ( !i ){
                         stop("Sample names are not matching 1:1 or misspecified in 'mtcount' and/or 'mtmeta'. Please read the documentation.")
                       }
                       # Check consistent gene names.
                       i <- setequal(mtgene[,1],mtcount[,1])
                       if ( !i ){
                         stop("Gene names are not matching 1:1 or misspecified in 'mtcount' and/or 'mtgene'. Please read the documentation.")
                       }

                       ### DUMP OUTPUT ###
                       sampOrd <- as.character(mtmeta$SampleID)

                       self$mtcount <- mtcount[,c(1,1+match(sampOrd,colnames(mtcount)[-1])),drop = F]
                       self$mtmeta  <- mtmeta[match(sampOrd,mtmeta$SampleID),,drop = F]
                       self$mtgene  <- mtgene[match(mtcount$GeneID,mtgene$GeneID),,drop = F]
                     },
                     ### WRITE THE PRINT FUNCTION --------------------------------------------------
                     print = function(...){

                       # Calculate row- and columnwise stats.
                       genesum <- self$mtcount[,-1] %>% as.matrix() %>%
                         base::rowSums(.) %>%
                         {c('Avg#reads' = mean(.,na.rm = T),
                            'min#reads' = min(.,na.rm = T),
                            'max#reads' = max(.,na.rm = T))} %>% round()
                       sampsum <- self$mtcount[,-1] %>% as.matrix() %>%
                         base::colSums(.) %>%
                         {c('Avg#reads' = mean(.,na.rm = T),
                            'min#reads' = min(.,na.rm = T),
                            'max#reads' = max(.,na.rm = T))} %>% round()
                       sums <- rbind(genesum,sampsum) %>% `rownames<-`(c("Genes:","Samples:"))

                       # Print
                       cat(
                         "mmt object with 3 elements.",
                         crayon::underline("\nSummary of count table:\n"))
                       print.table(sums,justify = "right")
                       cat(
                         crayon::underline("\nGenedata variables:"),ncol(self$mtgene),"\n",
                         paste(as.character(colnames(self$mtgene)),collapse = ", "),
                         crayon::underline("\nMetadata variables:"),ncol(self$mtmeta),"\n",
                         paste(as.character(colnames(self$mtmeta)),collapse = ", "))
                     }

                     ### WRITE THE SUBSET FUNCTION -------------------------------------------------
                     subset = function(mt,sub_genes = NULL,sub_samples = NULL,minreads = 0){

                       #For printing removed samples and OTUs
                       nsamplesbefore <- nrow(mt$meta) %>% as.numeric()
                       ngenesbefore   <- nrow(mt$count) %>% as.numeric()

                       #remove genes below minreads.
                       wh <- rowSums(mt$count[,-1]) >= minreads
                       mt$count <- mt$count[wh,,drop = F]
                       mt$gene  <- if(!is.null(mt$gene)) mt$gene <- mt$gene[wh,,drop = F]

                       # Subset samples.
                       if (!is.null(sub_samples)){
                         wh <- tryCatch(expr = {
                           mt$meta %>% filter(eval(parse(text = sub_samples)))
                         },error = function(e){
                           stop("The provided 'sub_samples' string is not meaningfull for the metadata.")
                         }) %>% {mt$meta$SampleID %in% .$SampleID}
                         mt$meta  <- mt$meta[wh,,drop = F]
                         mt$count <- mt$count[,c(TRUE,wh),drop = F]
                       }

                       # Subset genes.
                       if (!is.null(mt$gene)){
                         if(!is.null(sub_genes)){
                           wh <- tryCatch(expr = {
                             mt$gene %>% filter(eval(parse(text = sub_genes)))
                           },error = function(e){
                             stop("The provided 'sub_genes' string is not meaningfull for the gene data.")
                           }) %>% {mt$gene$GeneID %in% .$GeneID}
                           mt$count <- mt$count[wh,,drop = F]
                           mt$gene  <- mt$gene[wh,,drop = F]
                         }
                       } else {
                         stop("There is no gene data available for this mt object.")
                       }

                       #remove genes below minreads.
                       wh <- rowSums(mt$count[,-1]) >= minreads
                       mt$count <- mt$count[wh,,drop = F]
                       mt$gene  <- if(!is.null(mt$gene)) mt$gene <- mt$gene[wh,,drop = F]

                       # Report the removed samples/genes.
                       nsamplesafter <- nrow(mt$meta) %>% as.numeric()
                       ngenesafter   <- nrow(mt$count) %>% as.numeric()

                       message(paste(nsamplesbefore - nsamplesafter, "samples and",
                                     ngenesbefore - ngenesafter, "genes have been filtered \nBefore:",
                                     nsamplesbefore, "samples and", ngenesbefore, "genes\nAfter:",
                                     nsamplesafter, "samples and", ngenesafter, "genes"))
                     }


                   ))
