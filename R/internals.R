#' x without y
#'
"%w/o%" <- function(x, y) x[!x %in% y]

#' Prints mmt object summary (internal function)
#'
#' @param obj (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#'
#' @importFrom crayon underline
#' @method print mmt
#' @export
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
print.mmt <- function(obj){

  # Calculate row- and columnwise stats.
  genesum <- obj$mtdata[,-1] %>% as.matrix() %>%
    base::rowSums(.) %>%
    {c('Avg#reads'    = mean(.,na.rm = T),
       'Median#reads' = median(.,na.rm = T),
       'min#reads'    = min(.,na.rm = T),
       'max#reads'    = max(.,na.rm = T))} %>% round()
  sampsum <- obj$mtdata[,-1] %>% as.matrix() %>%
    base::colSums(.) %>%
    {c('Avg#reads'    = mean(.,na.rm = T),
       'Median#reads' = median(.,na.rm = T),
       'min#reads'    = min(.,na.rm = T),
       'max#reads'    = max(.,na.rm = T))} %>% round()
  sums <- rbind(genesum,sampsum) %>% `rownames<-`(c("Genes:","Samples:"))

  # Print
  cat(class(obj),"object with",nrow(obj$mtdata),"genes and",ncol(obj$mtdata)-1,"samples, consisting of 3 elements.",
    crayon::underline("\nSummary of data table:\n"))
  print.table(sums,justify = "right")
  cat(
    crayon::underline("\nGenedata variables:"),ncol(obj$mtgene),"\n",
    paste(as.character(colnames(obj$mtgene)),collapse = ", "),
    crayon::underline("\nMetadata variables:"),ncol(obj$mtmeta),"\n",
    paste(as.character(colnames(obj$mtmeta)),collapse = ", "),"\n")

  if (!is.null(attributes(obj)$normalised)){
    cat(crayon::underline("Normalised:\n"),
        paste0("Data is normalised by '",attributes(obj)$normalised,"'. See documentation for details.\n"))
  }
  if (!is.null(attributes(obj)$batch)){
    cat(crayon::underline("Batch-correction:\n"),
        paste0("Data is batch-corrected by variable '",attributes(obj)$batch$Var,"', using '",attributes(obj)$batch$method,"'. See documentation for details.\n"))
  }
}

#' Gather the mmt object into one dataframe, suitable for e.g. ggplot or custom analysis (internal function)
#'
#' @param mmt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param metavars (\emph{optional}) Columns from mtgene/mtmeta to keep. (\emph{default:} \code{NULL} (all columns))
#'
#' @importFrom dplyr inner_join one_of select everything
#' @importFrom data.table melt
#' @importFrom tidyr gather
#'
#' @return A data.table in long format. The gene count/expression is stored in the \code{Exprs} column.
#' @export
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
mt_gather <- function(mmt,metavars = c("SampleID","GeneID")){
  samps <- colnames(mmt$mtdata)[-1]

  # Setup the vectors for column selection (select all if metavars = NULL).
  if(is.null(metavars)) metavars <- c(colnames(mmt$mtgene),colnames(mmt$mtmeta))

  metaV <- metavars %w/o% c("SampleID",colnames(mmt$mtgene))
  geneV <- metavars %w/o% c("GeneID",colnames(mmt$mtmeta))
  if (!all(metaV %in% colnames(mmt$mtmeta)) | !all(geneV %in% colnames(mmt$mtgene)))
    stop("'metavars' must be valid column names from either mtgene or mtmeta.",call. = FALSE)

  # Join the mtgene and mtdata.
  gdat <- data.table:::merge.data.table(
    data.table:::subset.data.table(mmt$mtgene,select = c("GeneID",geneV)),
    mmt$mtdata,by = "GeneID")
  wh <- match(mmt$mtgene$GeneID,gdat$GeneID)
  setroworder(gdat,neworder = wh)
  gdat <- melt(gdat,
    id.vars       = colnames(gdat) %w/o% samps,
    measure.vars  = samps,
    value.name    = "Exprs",
    variable.name = "SampleID")
  gdat[,ID := .I]

  # Join with mtmeta.
  gdat <- data.table:::merge.data.table(
    data.table:::subset.data.table(mmt$mtmeta,select = c("SampleID",metaV)),
    gdat,by = "SampleID")
  setorder(gdat,ID)
  gdat[,ID := NULL]

  # Set correct column order.
  setcolorder(gdat,c("SampleID","GeneID","Exprs",metaV,geneV))
  return(gdat)
}
#' The row-version of setcolorder in data.table (internal function)
#'
#' @param x A \code{data.table}.
#' @param neworder Numeric vector of the row ordering.
#'
#' @import data.table
#'
#' @return A data.table
#' @export
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
setroworder <- function(x, neworder) {
    .Call(data.table:::Creorder, x, as.integer(neworder), PACKAGE = "data.table")
    invisible(x)
}
#' Computes a by-genome frequency matrix of any attribute in the mtgene data (internal function)
#'
#' @param mmt (\emph{required}) Data list as loaded with \code{\link{mt_load}}.
#' @param Genome (\emph{required}) Column name in mtgene grouping into genomes. NA and empty fields are ignored.
#' @param ID (\emph{required}) Column name in mtgene identifying the attribute. NA and empty fields are ignored.
#'
#' @import data.table
#'
#' @return A data.table
#' @export
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
mt_GenomeFrequencyMatrix <- function(mmt,Genome,ID){

  XX <- Genome
  YY <- ID

  out <- mmt$mtgene[,.(Freq = .N),by = mget(c(XX,YY))]
  setnames(out,old = 1:3,new = c("Genome","ID","Freq"))

  out <- out[!(Genome == "" | is.na(Genome) | ID == "" | is.na(Genome))]

  dcast(out,ID ~ Genome,fill = 0,value.var = "Freq")
}

#' performs a hypergeometric test on a set of genes (internal function)
#'
#' @param q_genes (\emph{required}) genes in query.
#' @param p_genes (\emph{required}) genes in pathway of interest.
#' @param u_genes (\emph{required}) genes in universe.
#' @param alternative (\emph{required}) test for enrichment ('greater'), depletion ('less') or both ('two.sided').
#'
#' @import data.table
#'
#' @return A named vector.
#' @export
#'
#' @author Thomas Yssing Michaelsen \email{tym@@bio.aau.dk}
test_pathway <- function(q_genes,p_genes,u_genes,alternative = "greater"){

  k  <- sum(q_genes %in% p_genes)
  n  <- length(q_genes)

  K  <- sum(u_genes %in% p_genes)
  N  <- length(u_genes)

  htest <- switch(alternative,
                  greater = {
                    Category:::.doHyperGInternal(
                      numW      = K,
                      numB      = N - K,
                      numDrawn  = n,
                      numWdrawn = k,
                      over      = T)
                  },
                  less = {
                    Category:::.doHyperGInternal(
                      numW      = K,
                      numB      = N - K,
                      numDrawn  = n,
                      numWdrawn = k,
                      over      = F)
                  },
                  two.sided = {
                    upr <- Category:::.doHyperGInternal(
                      numW      = K,
                      numB      = N - K,
                      numDrawn  = n,
                      numWdrawn = k,
                      over      = T)
                    lwr <- Category:::.doHyperGInternal(
                      numW      = K,
                      numB      = N - K,
                      numDrawn  = n,
                      numWdrawn = k,
                      over      = F)

                    list(p = {(2*min(c(upr$p,lwr$p))) %>% ifelse(. > 1,1,.)},odds = lwr$odds,expected = lwr$expected)
                  })

  c(
    unlist(htest),
    INquery    = k,
    Nquery     = n,
    INuniverse = K,
    Nuniverse  = N)
}

#' #' Generates an object of class 'diffexprs' (internal function)
#' #'
#' #'
#' #'
#' #'
#' diffexprs <- R6::R6Class("diffexprs",
#'   # PUBLIC #####################################################################
#'   public  = list(
#'     # Initialize with input data.
#'     initialize = function(dds){
#'       if (!class(dds) == "DESeqDataSet"){
#'         stop("You can only provide an object of class 'DESeqDataSet'")
#'       }
#'       if(length(as.character(dds@design)) > 2){
#'         stop("More than a single term in the model is not supported yet.")
#'       }
#'       private$..dds <- dds
#'     }
#'   ),
#'   # PRIVATE ####################################################################
#'   private = list(
#'     ..dds = NULL
#'   ),
#'   # ACTIVE #####################################################################
#'   active  = list(
#'     heatmapData = function(row_show = 10,row_label = NULL,row_order_by = NULL){
#'       #### INPUT CHECKING ###
#'       "%w/o%" <- function(x, y) x[!x %in% y]
#'       if(!is.null(row_label)){
#'         if(!(any(row_label %in% mcols(private$..dds)))){
#'           stop(paste0("Some/all of 'row_label' is not in the row-data."))
#'         } else {
#'           row_label <- c(colnames(varmeta)[1],row_label %w/o% colnames(varmeta)[1])
#'
#'
#'       wh <- which(colnames(mcols(private$..dds)) == "baseMean")
#'       if (wh == 1){
#'         id <- "ID"
#'       } else {
#'         id <- mcols(private$..dds) %>% {colnames(.)[1]}
#'       }
#'       # Get group and reference level.
#'       group <- private$..dds@design %>% {as.character(.)[2]}
#'       levs  <- {dds[[group]]} %>% levels()
#'
#'       # Run analysis for each level.
#'       res <- sapply((levs %w/o% levs[1]),function(i){
#'         results(private$..dds,contrast = c(group,levs[1],i))
#'       })
#'
#'       # Calculate logFC for each level in group.
#'       logFC <- lapply(names(res),function(i){
#'         x <- res[[i]] %>%
#'           data.frame() %>%
#'           select(log2FoldChange) %>%
#'           {eval(parse(text = paste0("dplyr::rename(.,'",i,"' = log2FoldChange)")))}
#'       }) %>% Reduce(f = function(x,y){merge(x,y,by = 0)},x = .) %>%
#'       {if (ncol(.) < 2)
#'         rownames_to_column(.,id)
#'         else
#'           rename(.,Row.names = id)}
#'
#'       # Order by variance and select top 'row_show'.
#'       if (ncol(logFC) < 3) {
#'         ord <- order(abs(logFC[,2]),decreasing = T)
#'       } else {
#'         ord <- order(apply(logFC[,2:ncol(logFC)],1,var),decreasing = T)
#'       }
#'       topFC <- logFC[ord,] %>%
#'       {if (nrow(.) <= row_show) .[1:nrow(.),] else .[1:row_show,]}
#'
#'       # Merge with mcols metadata.
#'       if(wh == 1){
#'
#'       } else {
#'         tmp <- mcols(dds)[,1:(wh-1)]
#'         logFC2 <- suppressWarnings(right_join(
#'           {tmp %>% select(row_label)},
#'           topFC,by = colnames(varmeta)[1])) %>%
#'           mutate(lab = paste(!!! rlang::syms(row_labels), sep = "; ")) %>%
#'           dplyr::select(-one_of(row_labels))
#'       }
#'
#'
#'
#'
#'     },
#'     ContrastData = function(){
#'
#'     },
#'     get_design = function(){
#'       private$..dds@design
#'     }
#' ))
#'
#'
#' de <- diffexprs$new(dds = dds)
#'
#' hest <- de$heatmapData
#' ko   <- de$.heatmap()


