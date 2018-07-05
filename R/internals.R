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
