#' mt_stackbar
#'
#' @title Make a stacked barplot of the sequencing statistics.
#'
#' @usage mt_stackbar(mmt,Raw,after_QC)
#'
#' @param mmt (\emph{required}) An object of class \code{mmt}.
#' @param size_axis (\emph{optional}) Specify size of primary and secondary axis. (\emph{Default: } \code{10})
#' @param size_text (\emph{optional}) Specify size of text in plot. (\emph{Default: } \code{3})
#' @param group_by  (\emph{optional}) A column name in the metadata to group samples by. (\emph{Default: } \code{NULL})
#' @param round_to  (\emph{optional}) Round the counts to thousands ("K") or milions ("M"). (\emph{Default: } \code{"M"})
#' @param detailed  (\emph{optional}) Output a table with the data used for the plot.
#' @return A ggplot
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom data.table data.table
#' @importFrom dplyr rename mutate select one_of left_join group_by starts_with distinct arrange
#' @importFrom tibble add_column
#' @importFrom tidyr gather
#' @importFrom plyr mapvalues
#'
#' @export
#'
#' @details The \code{\link{mt_stackbar}} function only works on data generated
#' from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline and loaded
#' using the \code{\link{mt_loadMetaT}} function. It is meant for data quality assesment
#' by giving a quick overview of the proportions of data removed for each pre-processing
#' step, if this is available in the metadata.
#'
#'
#' @examples
#'
#' \dontrun{
#' data("example_mmt")
#'
#' mt_stackbar(example_mmt,Raw = "Raw",after_QC = "after_QC",after_filter = "after_filter")
#' }
#'
#' @author Thomas Yssing Michaelsen \email{tym@bio.aau.dk}

mt_stackbar <- function(mmt,
                        group_by  = NULL,
                        round_to  = "M",
                        size_axis = 10,
                        size_text = 3,
                        detailed  = F){
  "%w/o%" <- function(x, y) x[!x %in% y]

  ##### CHECKING #####
  if(class(mmt) != "mmt") stop("Object is not of class 'mmt'")

  # Check data has not been normalised or batch-corrected.
  if(!is.null(attributes(mmt)$normalised)){
    stop("You cannot count raw # rRNA with normalised data.",call. = FALSE)
  }
  if(!is.null(attributes(mmt)$batch)){
    stop("You cannot count raw # rRNA with batch-corrected data.",call. = FALSE)
  }

  # Check the grouping variable.
  if(!is.null(group_by)){
    if (all(group_by %in% colnames(mmt$mtmeta))){
      if ("SampleID" %in% group_by){
        group_by <- c("SampleID",group_by %w/o% "SampleID")
      } else {
        group_by <- c("SampleID",group_by)
      }
      grp <- mmt$mtmeta[,.(
        SampleID = SampleID,
        Group    = Reduce(paste,.SD)),.SDcols = group_by]
    } else {
      stop("Your 'group_by' is not in metadata.",call. = FALSE)
    }
  } else {
    grp <- mmt$mtmeta[,.(SampleID = SampleID,Group = as.character(SampleID))]
  }

  # Check if seqstat data are available.
  wh <- c("Raw","after_QC")
  i  <- sapply(wh,match,table = colnames(mmt$mtmeta))
  if(any(is.na(i))){
    err <- names(i[which(is.na(i))]) %>% paste(collapse = ", ")
    stop(paste0("The following is missing from metadata: ",err),call. = FALSE)
  }

  # Accumulate rRNA and non-rRNA counts.
  tmp <- mt_gather(mmt,metavars = "ftype")[,.(
    SampleID = SampleID,GeneID = GeneID,Exprs = Exprs,
    ftype    = ifelse(ftype == "rRNA","rRNA","nonrRNA"))][,.(
    x        = sum(Exprs)),by = .(SampleID,ftype)] %>%
    dcast(SampleID ~ ftype,value.var = "x") %>%
    .[,.(
      SampleID = SampleID,rRNA = rRNA,nonrRNA = nonrRNA,
      mapped = nonrRNA+rRNA)] %>%
    setkey(SampleID)

  mmt$mtmeta <- tmp[mmt$mtmeta]

  # Check values are decreasing.
  dat <- mmt$mtmeta %>% select(Raw,after_QC,mapped)
  i   <- apply(dat,1,function(x){all(diff(as.numeric(x)) <= 0)})
  if (any(!i)){
    err <- mmt$mtmeta[i,1,drop = T] %>% as.character() %>% paste(collapse = ", ")
    stop(paste0("A decrease like Raw >= after_QC >= mapped was not found in: ",err),call. = FALSE)
  }

  ##### PREPARE DATA #####
  tab <- as.data.frame(mmt$mtmeta,stringsAsFactors = F) %>%
    select(Raw,after_QC,mapped,nonrRNA) %>%
    apply(.,1,function(x){
      abs(diff(x)) %>% c(.,Left = unname(x[1]) - sum(.)) %>% `names<-`(.,paste0("raw_",names(.)))
    }) %>%
    t() %>% as.data.frame() %>%
    {cbind(.,t(apply(.,1,function(x){`names<-`(prop.table(x)*100,paste0("pct_",colnames(.)))})))} %>%
    add_column(SampleID = mmt$mtmeta$SampleID,.before = 1) %>%
    left_join(.,grp,by = "SampleID") %>%
    {left_join(
      select(.,SampleID,Group,dplyr::starts_with("pct_raw_")) %>% gather(key = RNA,value = Freq,-SampleID,-Group) %>% mutate(RNA = gsub("pct_raw_","",RNA)),
      select(.,SampleID,Group,dplyr::starts_with("raw_")) %>% gather(key = RNA,value = Cnt,-SampleID,-Group) %>% mutate(RNA = gsub("raw_","",RNA)),
      by = c("SampleID","Group","RNA"))} %>%
    mutate(RNA = factor(RNA,levels = unique(RNA))) %>%
    mutate(Group = factor(Group,levels = unique(Group))) %>%
    dplyr::group_by(Group)

  # Round the counts.
  rnd <- switch(round_to,
                K = 10^3,
                M = 10^6)
  if(is.null(rnd)) stop("Please specify 'round_to' correctly. See documentation",call. = FALSE)

  tab <- tab %>%
    mutate(Total = paste(round(sum(Cnt)/rnd,2),round_to)) %>%
    mutate(Cnt.raw = Cnt) %>%
    mutate(Cnt     = round(Cnt/rnd,2)) %>%
    mutate(lab     = 100) %>%
    mutate(RNA     = plyr::mapvalues(RNA,
                                     from = c("after_QC","mapped","nonrRNA","Left"),
                                     to   = c("Removed by QC","unmapped","rRNA","non-rRNA")))

  # Top axis labels.
  llab <- tab %>% dplyr::distinct(Group,Total) %>% arrange(Group)

  ##### PLOTTING #####
  p <- ggplot(tab, aes(x = as.numeric(Group), y = Freq, fill = RNA,label = Cnt,text = Total,text_pos = lab)) +
    geom_bar(stat = "identity") +
    geom_text(size = size_text, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(breaks=seq(0,100,10)) +
    scale_x_continuous(
      breaks   = 1:nrow(llab),
      labels   = llab$Group,
      sec.axis = sec_axis(~.,breaks = 1:nrow(llab),labels = llab$Total)) +
    scale_fill_manual(values = c("red","grey","yellow","green")) +
    theme(
      axis.text.x      = element_text(size = size_axis, color = "black", angle = 90, hjust = 1, vjust = 0.4),
      axis.title.x     = element_blank(),
      axis.line        = element_blank(),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      legend.title     = element_blank(),
      axis.text.x.top  = element_text(angle = 45,face = "bold",hjust = .5,size = size_axis*.8),
      axis.ticks.x     = element_blank()) +
    ylab("Composition (%)")

  if(!detailed){
    return(p)
  } else {
    # Output table.
    out <- mmt$mtmeta %>% select(SampleID,Raw,after_QC,mapped,rRNA)
    return(list(plot = p,table = out))
  }
}
