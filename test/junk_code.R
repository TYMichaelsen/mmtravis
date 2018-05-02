## Junk code that may be handy later.




mt_stats <- function(mt,group_by = NULL){
  "%w/o%" <- function(x, y) x[!x %in% y]

  # Check the grouping variable.
  if(!is.null(group_by)){
    if (all(group_by %in% colnames(mt$meta))){
      if ("SampleID" %in% group_by){
        group_by <- c("SampleID",group_by %w/o% "SampleID")
      } else {
        group_by <- c("SampleID",group_by)
      }
      grp <- cbind(mt$meta[1],Group = apply(mt$meta[,group_by,drop = F], 1, paste, collapse = " ")) %>%
        mutate(Group = as.character(Group))
    } else {
      stop("Your 'group_by' is not in metadata.",call. = FALSE)
    }
  } else {
    grp <- cbind(mt$meta[1],Group = mt$meta[,1]) %>%
      mutate(Group = as.character(Group))
  }

  if(!is.null(mt$stat)){
    # Create output table.
    tab <- mt$stat %>%
      select(-SampleID) %>%
      apply(.,1,function(x){
        abs(diff(x)) %>% c(.,Left = unname(x[1]) - sum(.)) %>% `names<-`(.,paste0("raw_",names(.)))
      }) %>%
      t() %>% as.data.frame() %>%
      {cbind(.,t(apply(.,1,function(x){`names<-`(prop.table(x)*100,paste0("pct_",colnames(.)))})))} %>%
      add_column(SampleID = mt$stat$SampleID,.before = 1) %>%
      left_join(.,grp,by = "SampleID") %>%
      {left_join(
        select(.,SampleID,Group,dplyr::starts_with("pct_raw_")) %>% gather(key = RNA,value = Freq,-SampleID,-Group) %>% mutate(RNA = gsub("pct_raw_","",RNA)),
        select(.,SampleID,Group,dplyr::starts_with("raw_")) %>% gather(key = RNA,value = Cnt,-SampleID,-Group) %>% mutate(RNA = gsub("raw_","",RNA)),
        by = c("SampleID","Group","RNA"))} %>%
      mutate(RNA = factor(RNA,levels = unique(RNA))) %>%
      mutate(Group = factor(Group,levels = unique(Group))) %>%
      dplyr::group_by(Group) %>%
      dplyr::mutate(Total = paste(round(sum(Cnt)/10^6,2),"mio")) %>%
      dplyr::mutate(Cnt   = round(Cnt/10^6,2)) %>%
      mutate(lab   = 100) %>%
      mutate(RNA   = mapvalues(RNA,
                               from = c("QC","filtered","Left"),
                               to   = c("Removed by QC","Removed rRNA","Used for analysis")))

    # Prep stuff for plotting.
    llab <- tab %>% distinct(Group,Total) %>% arrange(Group)

    # Plot the output.
    ggplot(tab, aes(x = as.numeric(Group), y = Freq, fill = RNA,label = Cnt,text = Total,text_pos = lab)) +
      geom_bar(stat = "identity") +
      geom_text(size = 3, position = position_stack(vjust = 0.5)) +
      scale_y_continuous(breaks=seq(0,100,10)) +
      scale_x_continuous(
        breaks   = 1:nrow(llab),
        labels   = llab$Group,
        sec.axis = sec_axis(~.,breaks = 1:nrow(llab),labels = llab$Total)) +
      theme(
        axis.text.x      = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.4),
        axis.title.x     = element_blank(),
        axis.line        = element_blank(),
        panel.background = element_blank(),
        panel.border     = element_blank(),
        legend.title     = element_blank(),
        axis.text.x.top  = element_text(angle = 45,face = "bold",hjust = .5,size = 8),
        axis.ticks.x     = element_blank()) +
      ylab("Composition (%)")

  } else {
    stop("You do not have a 'mtstat' object assigned to your data.",call. = FALSE)
  }
}
