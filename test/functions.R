################################################################################
# Functions for transcriptome analysis.
################################################################################

################################################################################
trans_dumpDE <- function(DEobj,DEres,prefix = "",path){
  # Build table.
  tab <- merge({
   DEres %>%
      as.data.frame() %>%
      rownames_to_column("ID") %>%
      {add_column(.,FoldChange = 2^(.$log2FoldChange),.after = "log2FoldChange")}
  },
  {
    counts(DEobj,normalized = T) %>%
      as.data.frame() %>%
      rownames_to_column("ID")
  },by = "ID") %>%
    merge({
    mcols(DEobj) %>%
      as.data.frame() %>%
      {.[,-c(which(names(.) == "baseMean"):ncol(.))]}
  },.,by = "ID",all.y = T) %>%
    subset(.,padj < 0.05) %>%
    arrange(.,desc(FoldChange))

  # Extract group names for file.
  nm <- DEres@elementMetadata$description %>%
    grep(pattern = "log2 fold change",value = T) %>%
    {strsplit(x = .,split = ": ")[[1]][2]} %>%
    {strsplit(x = .,split = " ")[[1]][c(1,2,4)]} %>%
    {paste0(.[2],"_vs_",.[3])}

  # Write file.
  write.table(
    x         = tab,
    file      = paste0(path,"/",prefix,nm,".txt"),
    row.names = F,
    quote     = F,
    sep       = "\t")
}

################################################################################

trans_MAplot <- function(res,rename_contrast = NULL){
  # Fetch the contrast names from object.
  if(is.null(rename_contrast)){
    ctr <- res@elementMetadata$description %>%
      grep(pattern = "log2 fold change",value = T) %>%
      {strsplit(x = .,split = ": ")[[1]][2]} %>%
      {strsplit(x = .,split = " ")[[1]][c(1,2,4)]} %>%
      {c(paste0(.[1]," = ",.[2]),paste0(.[1]," = ",.[3]))}
  } else {
    ctr <- rename_contrast
  }

  # Convert to data.frame for plotting.
  res <- res %>%
  {cbind.data.frame(ID = rownames(.),.)} %>%
    mutate(Significant = ifelse(padj < 0.01, "Yes", "No"),
           ER          = ifelse(log2FoldChange > 0,ctr[1],ctr[2]))

  # Count the number of differentially expressed.
  labs <- subset(res,select = c(Significant,ER)) %>%
    table %>% .[2,ctr] %>% {paste(names(.),.,sep = "\nn = ")}

  ## plot.
  ggplot(res, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
    geom_hline(yintercept = 0, color = "darkred", lty = 2) +
    geom_point(size = 1) +
    scale_color_manual(na.value = "black", values = c("black", "red")) +
    theme_classic() +
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          legend.position = "none") +
    annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, label = labs[1], size = 3) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -.2, label = labs[2], size = 3) +
    ylab("Fold Change (log2)") +
    xlab("Expression (base mean)")
}

################################################################################
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.text <- function(x, y,txt,...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))

  text(0.5, 0.5, txt,...)
}

################################################################################
glOverviewPlot <- function(res,DEobj,group,...){

  # glXYPlot cannot handle missing values, remove them.
  wh <- which(!is.na(res$padj))
  res_sub   <- res[wh,]
  DEobj_sub <- DEobj[wh,]
  sig_sub   <- as.numeric(res_sub$padj < 0.05)

  glXYPlot(
    x         = log10(res_sub$baseMean),
    y         = -log10(res_sub$padj),
    xlab      = "logMean",
    ylab      = "padj",
    side.xlab = group,
    status    = sig_sub,
    counts    = counts(DEobj_sub,normalized = T),
    groups    = DEobj_sub$Time,
    samples   = colnames(DEobj_sub),
    anno      = mcols(DEobj_sub)[,c("ID","contig","ftype","gene","product")],
    side.main = "ID",...)
}

################################################################################
concat <- function(max.length = 50,x,y){
  mapply(function(x,y){
    nx <- nchar(x %>% as.character)
    ny <- nchar(y %>% as.character)
    if (sum(nx,ny) > max.length){
      paste0(substr(x,1,(max.length-5)-ny),"... ;",y)
    } else {
      paste(x,y,sep = "; ")
    }
  },x = x,y = y)
}

################################################################################
# ```{r,echo=F,warning=F,eval=F}
#
# sig <- as.numeric(res$padj<0.05)
# glMDPlot(res,
#          status    = sig,
#          counts    = counts(DEobj,normalized = T),
#          groups    = DEobj$Time,
#          transform = FALSE,
#          samples   = colnames(DEobj),
#          anno      = mcols(DEobj)[,c("ID","contig","ftype","gene","product")],
#          side.main = "ID")
#
# ```

################################################################################
trans_ordinate <- function(data,
  filter_genes = 0,
  type = "pca",
  distmeasure = "euclidian",
  transform = "none",
  constrain = NULL,
  x_axis = 1,
  y_axis = 2,
  sample_color_by = NULL,
  sample_color_order = NULL,
  sample_shape_by = NULL,
  sample_colorframe = FALSE,
  sample_colorframe_label = NULL,
  sample_label_by = NULL,
  sample_label_size = 4,
  sample_label_segment_color = "black",
  sample_point_size = 2,
  sample_trajectory = NULL,
  sample_trajectory_group = sample_trajectory,
  sample_plotly = NULL,
  genes_plot = FALSE,
  genes_nlabels = 0,
  genes_label_by = "ID",
  genes_label_size = 3,
  genes_label_color = "grey10",
  genes_rescale = FALSE,
  genes_point_size = 2,
  genes_shape = 20,
  genes_plotly = FALSE,
  envfit_factor = NULL,
  envfit_numeric = NULL,
  envfit_signif_level = 0.005,
  envfit_textsize = 3,
  envfit_textcolor = "darkred",
  envfit_numeric_arrows_scale = 1,
  envfit_arrowcolor = "darkred",
  envfit_show = TRUE,
  repel_labels = TRUE,
  opacity = 0.8,
  detailed_output = FALSE,
  ...) {

  ##### Sanity check of options  #####
  if(genes_plotly == TRUE & !is.null(sample_plotly)){
    stop("You can not use plotly for both genes and samples in the same plot.", call. = FALSE)
  }
  if(genes_plotly == TRUE | !is.null(sample_plotly)){
    #message("geom_text_repel is not supported by plotly yet.")
    repel_labels <- FALSE
  }

  #Impossible to do ordination with 1 or 2 samples
  if(length(unique(data$metadata[,1])) <= 2)
    stop("Ordination cannot be performed on 2 or fewer samples (the number of resulting axes will always be n-1, where n is the number of samples).", call. = FALSE)

  if(is.null(sample_color_by) & !is.logical(sample_colorframe) & !is.null(sample_colorframe)) {
    sample_color_by <- sample_colorframe
  }

  ##### Filter #####
  if(filter_genes != 0) {
    #First transform to percentages
    abund_pct <- data$exprs %>% as.data.frame()
    abund_pct[,which(colSums(abund_pct) != 0)] <- as.data.frame(apply(abund_pct[,which(colSums(abund_pct) != 0), drop = FALSE], 2, function(x) x/sum(x)*100))
    rownames(abund_pct) <- rownames(data$exprs) #keep rownames

    #Then filter low abundant genes where ALL samples have below the threshold set with filter_genes in percent
    abund_subset <- abund_pct[!apply(abund_pct, 1, function(row) all(row <= filter_genes)),,drop = FALSE] #remove low abundant genes
    data$exprs <- data$exprs[which(rownames(data$exprs) %in% rownames(abund_subset)),,drop = FALSE]
    rownames(data$genedata) <- data$genedata$ID
    data$genedata <- data$genedata[which(rownames(data$genedata) %in% rownames(abund_subset)),,drop = FALSE] #same with taxonomy
  }

  #to fix user argument characters, so fx PCoA/PCOA/pcoa are all valid
  type <- tolower(type)

  ##### Data transformation with decostand()  #####
  if(!transform == "none" & transform != "sqrt") {
    transform <- tolower(transform)
    data$exprs <- t(vegan::decostand(t(data$exprs), method = transform))
  } else if (tolower(transform) == "sqrt") {
    data$exprs <- t(sqrt(t(data$exprs)))
  }

  ##### Inputmatrix AFTER transformation  #####
  if (any(type == c("nmds", "mmds", "pcoa", "dca"))) {
    if(!type == "nmds" & (genes_plot == TRUE | genes_plotly == TRUE)) {
      stop("No speciesscores available with mMDS/PCoA, DCA.", call. = FALSE)
    }
    if (!distmeasure == "none") {
      #Calculate distance matrix with vegdist()
      distmeasure <- tolower(distmeasure)
      if(distmeasure == "jsd") {
        #This is based on http://enterotype.embl.de/enterotypes.html
        #Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
        #Unfortunately this code is SLOOOOOOOOW
        dist.JSD <- function(inMatrix, pseudocount=0.000001) {
          KLD <- function(x,y) sum(x *log(x/y))
          JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
          matrixColSize <- length(colnames(inMatrix))
          matrixRowSize <- length(rownames(inMatrix))
          colnames <- colnames(inMatrix)
          resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

          inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

          for(i in 1:matrixColSize) {
            for(j in 1:matrixColSize) {
              resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                                     as.vector(inMatrix[,j]))
            }
          }
          colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
          as.dist(resultsMatrix)->resultsMatrix
          attr(resultsMatrix, "method") <- "dist"
          return(resultsMatrix)
        }
        message("Calculating Jensen-Shannon Divergence (JSD) distances... ")
        inputmatrix <- dist.JSD(data$exprs)
        message("Done.")
      } else if(any(distmeasure == c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"))) {
        message("Calculating distance matrix... ")
        inputmatrix <- vegan::vegdist(t(data$exprs), method = distmeasure)
        message("Done.")
      }
    } else if (distmeasure == "none") {
      warning("No distance measure selected, using raw data. If this is not deliberate, please provide one with the argument: distmeasure.", call. = FALSE)
      inputmatrix <- t(data$exprs)
    }

    if (transform != "none" & distmeasure != "none") {
      warning("Using both transformation AND a distance measure is not recommended for distance-based ordination (nMDS/PCoA/DCA). If this is not deliberate, consider transform = \"none\".", call. = FALSE)
    }
  } else if(any(type == c("pca", "rda", "ca", "cca"))) {
    inputmatrix <- t(data$exprs)
  }

  ##### Perform ordination  #####
  #Generate data depending on the chosen ordination type
  if(type == "pca") {
    #make the model
    model <- vegan::rda(inputmatrix)

    #axis (and data column) names
    x_axis_name <- paste0("PC", x_axis)
    y_axis_name <- paste0("PC", y_axis)

    #Calculate the amount of inertia explained by each axis
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)

    #Calculate species- and site scores
    samplescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    genescores   <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "rda") {
    if(is.null(constrain))
      stop("Argument constrain must be provided when performing constrained/canonical analysis.", call. = FALSE)
    #make the model
    codestring <- paste0("rda(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))

    #axes depend on the results
    x_axis_name <- paste0("RDA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "PC1"

      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("RDA", y_axis)

      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }

    #Calculate species- and site scores
    samplescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    genescores   <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "nmds") {
    #make the model
    if(ncol(data$exprs) > 100) {
      message("Performing non-Metric Multidimensional Scaling on more than 100 samples, this may take some time ... ")
    }
    model <- vegan::metaMDS(inputmatrix, trace = FALSE, ...)
    if(ncol(data$exprs) > 100) {
      message("Done.")
    }

    #axis (and data column) names
    x_axis_name <- paste0("NMDS", x_axis)
    y_axis_name <- paste0("NMDS", y_axis)

    #Calculate species- and site scores
    #Speciesscores may not be available with MDS
    samplescores <- vegan::scores(model, display = "sites")
    if(!length(model$species) > 1) {
      genescores <- NULL
      if(sample_plot == TRUE | genes_plotly == TRUE) {
        stop("Genescores are not available.", call. = FALSE)
      }
    } else {
      samplescores <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
    }
  } else if(type == "mmds" | type == "pcoa") {
    #make the model
    model <- ape::pcoa(inputmatrix, ...)

    #axis (and data column) names
    x_axis_name <- paste0("PCo", x_axis)
    y_axis_name <- paste0("PCo", y_axis)

    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$values$Relative_eig * 100, 1)
    names(totalvar) <- c(paste0("PCo", seq(1:length(totalvar))))

    #Calculate species- and site scores
    #Speciesscores are not available with pcoa
    sitescores <- as.data.frame(model$vectors)
    colnames(sitescores) <- c(paste0("PCo", seq(1:length(samplescores))))
    genescores <- NULL
  } else if(type == "ca") {
    #make the model
    model <- vegan::cca(inputmatrix, ...)

    #axis (and data column) names
    x_axis_name <- paste0("CA", x_axis)
    y_axis_name <- paste0("CA", y_axis)

    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)

    #Calculate species- and site scores
    samplescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    genescores   <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "cca") {
    if(is.null(constrain))
      stop("Argument constrain must be provided when performing constrained/canonical analysis.", call. = FALSE)
    #make the model
    codestring <- paste0("cca(inputmatrix~", paste(constrain, collapse = "+"), ", data$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))

    #axes depend on the results
    x_axis_name <- paste0("CCA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "CA1"

      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("CCA", y_axis)

      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }

    #Calculate species- and site scores
    samplescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    genescores   <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "dca") {
    #make the model
    model <- vegan::decorana(inputmatrix, ...)

    #axis (and data column) names
    x_axis_name <- paste0("DCA", x_axis)
    y_axis_name <- paste0("DCA", y_axis)

    #Calculate the percentage of eigenvalues explained by the axes
    #totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)

    #Calculate species- and site scores
    samplescores <- vegan::scores(model, display = "sites", choices = c(x_axis, y_axis))
    genescores   <- vegan::scores(model, display = "species", choices = c(x_axis, y_axis))
  }

  ##### Data for ggplot  #####
  dsample <- cbind.data.frame(data$metadata, samplescores)

  if (!is.null(sample_color_order)) {
    dsample[, sample_color_by] <- factor(dsample[, sample_color_by], levels = sample_color_order)
  }

  if(length(genescores) > 1) {
    dgene <- merge(data.frame(genescores, ID = rownames(genescores)), data$genedata, by.x = "ID")
    dgene$dist <- dgene[, x_axis_name]^2 + dgene[, y_axis_name]^2
    dgene <- dplyr::arrange(dgene, desc(dist))
    rownames(dgene) <- dgene$ID
    if (genes_rescale == TRUE) {
      maxx <- max(abs(dsample[, x_axis_name]))/max(abs(dgene[,x_axis_name]))
      dspecies[, x_axis_name] <- dsample[, x_axis_name] * maxx * 0.8
      maxy <- max(abs(dsample[, y_axis_name]))/max(abs(dgene[,y_axis_name]))
      dgene[, y_axis_name] <- dgene[, y_axis_name] * maxy * 0.8
    }
  } else {
    dsample = NULL
  }

  ##### Base plot object #####
  plot <- ggplot(dsample,
                 aes_string(x = x_axis_name,
                            y = y_axis_name,
                            color = sample_color_by,
                            shape = sample_shape_by))

  ##### Colorframe  #####
  if(!sample_colorframe == FALSE) {
    if(is.null(sample_color_by) & sample_colorframe == TRUE)
      stop("Applying a colorframe to the sample points requires a variable in the metadata to be provided by the argument sample_colorframe, or by sample_color_by if the former is only TRUE.", call. = FALSE)
    if(sample_colorframe == TRUE) {
      splitData <- base::split(plot$data, plot$data[, sample_color_by]) %>%
        lapply(function(df) {
          df[chull(df[, x_axis_name], df[, y_axis_name]), ]
        })
      hulls <- do.call(rbind, splitData)
      plot <- plot + geom_polygon(data = hulls, aes_string(fill = sample_color_by, group = sample_color_by), alpha = 0.2*opacity)
    } else if (!is.logical(sample_colorframe) & !is.null(sample_colorframe)) {
      plot$data$colorframeGroup <- paste(plot$data[,sample_color_by], plot$data[,sample_colorframe]) %>% as.factor()
      splitData <- base::split(plot$data, plot$data$colorframeGroup) %>%
        lapply(function(df) {
          df[chull(df[, x_axis_name], df[, y_axis_name]), ]
        })
      hulls <- do.call(rbind, splitData)
      plot <- plot + geom_polygon(data = hulls, aes_string(fill = sample_color_by, group = "colorframeGroup"), alpha = 0.2*opacity)
    }
  }

  ##### Plot sample points  #####
  if (!is.null(sample_plotly)){
    if(length(sample_plotly) > 1){
      data_plotly <- apply(data$metadata[,sample_plotly], 1, paste, collapse = "<br>")
    } else if(sample_plotly == "all" | sample_plotly == TRUE){
      data_plotly <- apply(data$metadata[,], 1, paste, collapse = "<br>")
    } else{
      data_plotly <- paste0(sample_plotly,": ",data$metadata[,sample_plotly])
    }
    plot <- plot +
      suppressWarnings(geom_point(size = 2, alpha = opacity,
                                  aes(text = data_plotly))) + #HER
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  } else{
    plot <- plot +
      geom_point(size = sample_point_size, alpha = opacity) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "black", size = 0.5))
  }

  #Only eigenvalue-based ordination methods can be displayed with % on axes
  if(type == "pca" | type == "ca" | type == "pcoa" | type == "mmds") {
    plot <- plot +
      xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
      ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
  } else if(type == "rda" | type == "cca") {
    if(model$CCA$rank > 1) {
      plot <- plot +
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "% / ", constrainedvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    } else if(model$CCA$rank <= 1) {
      plot <- plot +
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    }
  } else if (type == "nmds") {
    plot <- plot +
      annotate("text", size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste0("Stress value = ", round(model$stress, 3)))
  }

  ##### Plot gene points  #####
  if (genes_plot == TRUE) {
    # if(genes_plotly == T){
    #   data_plotly <- paste("Kingdom: ", data$tax[,1],"<br>",
    #                        "Phylum: ", data$tax[,2],"<br>",
    #                        "Class: ", data$tax[,3],"<br>",
    #                        "Order: ", data$tax[,4],"<br>",
    #                        "Family: ", data$tax[,5],"<br>",
    #                        "Genus: ", data$tax[,6],"<br>",
    #                        "Species: ", data$tax[,7],"<br>",
    #                        "OTU: ", data$tax[,8],sep = "")
    #   plot <- plot +
    #     geom_point(data = dspecies,
    #                color = "darkgrey",
    #                shape = species_shape,
    #                size = species_point_size-1,
    #                alpha = opacity,
    #                aes(text = data_plotly))
    # } else{
    #   plot <- plot +
    #     geom_point(data = dspecies,
    #                color = "darkgrey",
    #                shape = species_shape,
    #                size = species_point_size,
    #                alpha = opacity)
    # }
  }

  ##### Plot text labels  #####
  if (!is.null(sample_colorframe_label)) {
    temp <- data.frame(group = dsample[, sample_colorframe_label],
                       x = dsample[, x_axis_name],
                       y = dsample[, y_axis_name]) %>%
      group_by(group) %>%
      summarise(cx = mean(x), cy = mean(y)) %>%
      as.data.frame()
    temp2 <- merge(dsample, temp,
                   by.x = sample_colorframe_label,
                   by.y = "group")
    temp3 <- temp2[!duplicated(temp2[, sample_colorframe_label]), ]
    if (repel_labels == T){
      plot <- plot + ggrepel::geom_text_repel(data = temp3,
                                              aes_string(x = "cx",
                                                         y = "cy",
                                                         label = sample_colorframe_label),
                                              size = 3,
                                              color = "black",
                                              fontface = 2
      )
    } else {
      plot <- plot + geom_text(data = temp3,
                               aes_string(x = "cx",
                                          y = "cy",
                                          label = sample_colorframe_label),
                               size = 3,
                               color = "black",
                               fontface = 2)
    }
  }

  ##### Sample_trajectory  #####
  if (!is.null(sample_trajectory)) {
    traj <- dsample[order(dsample[, sample_trajectory]), ]
    plot <- plot + geom_path(data = traj, aes_string(group = sample_trajectory_group))
  }

  ##### Sample point labels  #####
  if(!is.null(sample_label_by)) {

    if (repel_labels == T){
      plot <- plot +
        ggrepel::geom_text_repel(
          aes_string(label = sample_label_by),
          size          = sample_label_size,
          color         = "grey40",
          segment.color = sample_label_segment_color)}
    else{
      plot <- plot +
        geom_text(
          aes_string(label = sample_label_by),
          size          = sample_label_size,
          color         = "grey40",
          segment.color = sample_label_segment_color)}
  }

  ##### Plot gene labels  #####
  if (genes_nlabels > 0) {
    if (repel_labels == T){
      plot <- plot +
        ggrepel::geom_text_repel(data = dgene[1:genes_nlabels,],
          aes_string(x     = x_axis_name,
                     y     = y_axis_name,
                     label = genes_label_by),
          colour   = genes_label_color,
          size     = genes_label_size,
          fontface = 4,
          inherit.aes = FALSE)}
    else{
      plot <- plot +
        geom_text(data = dgene[1:genes_nlabels,],
          aes_string(x     = x_axis_name,
                     y     = y_axis_name,
                     label = genes_label_by),
          colour      = genes_label_color,
          size        = genes_label_size,
          fontface    = 4,
          inherit.aes = FALSE)}
  }

  ##### Categorical fitting  #####
  if(!is.null(envfit_factor)) {
    evf_factor_model <- envfit(dsample[,c(x_axis_name, y_axis_name)],
                               data$metadata[,envfit_factor, drop = FALSE],
                               permutations = 999
    )
    evf_factor_data <- data.frame(Name = rownames(evf_factor_model$factors$centroids),
                                  Variable = evf_factor_model$factors$var.id,
                                  evf_factor_model$factors$centroids,
                                  pval = evf_factor_model$factors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_factor_data) > 0 & envfit_show == TRUE) {
      if (repel_labels == T){plot <- plot + ggrepel::geom_text_repel(data = evf_factor_data,aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_textcolor, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")}
      else{plot <- plot + geom_text(data = evf_factor_data,aes_string(x = x_axis_name, y = y_axis_name, label = "Name"), colour = envfit_textcolor, inherit.aes = FALSE, size = envfit_textsize, fontface = "bold")}
    }
    if (nrow(evf_factor_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.", call. = FALSE)
    }
  } else {
    evf_factor_model <- NULL
  }

  ##### Numerical fitting  #####
  if (!is.null(envfit_numeric)) {
    evf_numeric_model <- envfit(dsample[,c(x_axis_name, y_axis_name)],
                                data$metadata[,envfit_numeric, drop = FALSE],
                                permutations = 999
    )
    evf_numeric_data <- data.frame(Name = rownames(evf_numeric_model$vectors$arrows),
                                   evf_numeric_model$vectors$arrows * sqrt(evf_numeric_model$vectors$r) * envfit_numeric_arrows_scale,
                                   pval = evf_numeric_model$vectors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_numeric_data) > 0 & envfit_show == TRUE) {
      plot <- plot + geom_segment(data = evf_numeric_data,
                                  aes_string(x = 0,
                                             xend = x_axis_name,
                                             y = 0,
                                             yend = y_axis_name
                                  ),
                                  arrow = arrow(length = unit(3, "mm")),
                                  colour = envfit_arrowcolor,
                                  size = 1,
                                  inherit.aes = FALSE) +
        geom_text(data = evf_numeric_data,
                  aes_string(x = x_axis_name,
                             y = y_axis_name,
                             label = "Name"),
                  colour = envfit_textcolor,
                  inherit.aes = FALSE,
                  size = envfit_textsize,
                  hjust = 1.2,
                  vjust = 1.2,
                  fontface = "bold")
    }
    if (nrow(evf_numeric_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.", call. = FALSE)
    }
  } else {
    evf_numeric_model <- NULL
  }

  ##### Return  #####
  #return plot or additional details
  if(!is.null(sample_plotly)){
    plotly::ggplotly(plot, tooltip = "text")
  }
  else if(genes_plotly == T){
    plotly::ggplotly(plot, tooltip = "text")
  }
  else if(!detailed_output){
    return(plot)
  }
}

################################################################################
trans_heatmap <- function(data,
group_by = NULL,
facet_by = NULL,
gene_show = 10,
sort_rows_by   = NULL,
sort_rows_func = var,
gene_aggregate = "ID",
gene_aggregate_func = mean,
sample_aggregate = NULL,
sample_aggregate_func = mean,
plot_colorscale = "log10",
color_vector = NULL,
plot_values = TRUE,
plot_values_size = 4,
round = 1,
sort_by = NULL,
rel_widths = c(0.75, 0.25)
) {

  ## Extract the data into separate objects for readability
  exprs    <- data[["exprs"]]
  genedata <- data[["genedata"]]
  metadata <- data[["metadata"]]

  ## Coerce the group_by and facet_by variables to factor to always be considered categorical. Fx Year is automatically loaded as numeric by R, but it should be considered categorical.
  ## Grouping a heatmap by a continuous variable doesn't make sense
  if(!is.null(group_by)) {
    metadata[group_by] <- lapply(metadata[group_by], factor)
  }

  if(!is.null(facet_by)) {
    if(is.null(group_by)) {
      group_by <- facet_by
    }
    metadata[facet_by] <- lapply(metadata[facet_by], factor)
  }

  # select a specific variable in the genedata to aggregate into.
  exprs3 <- cbind.data.frame(Display = genedata[[gene_aggregate]], exprs) %>%
    tidyr::gather(key = Sample, value = Exprs, -Display) %>% as.data.table()

  exprs3 <- exprs3[,
    aggr_gene:=gene_aggregate_func(Exprs),
    by   =list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    as.data.frame()

  if(!is.null(facet_by)){
    ogroup <- group_by
    group_by <- c(group_by, facet_by)
  }

  # Merge grouping variables with data.
  suppressWarnings(
    if (!is.null(group_by)){
      if (length(group_by) > 1){
        oldGroup <- unique(cbind.data.frame(
          Group  = apply(metadata[,group_by], 1, paste, collapse = " "),
          metadata[,group_by]))
        exprs3 <- merge(exprs3,
                        data.frame(Sample = metadata[,1],oldGroup),
                        by = "Sample")
      } else{
        grp <- data.frame(Sample = metadata[,1], Group = metadata[,group_by])
        exprs3$Group <- grp$Group[match(exprs3$Sample, grp$Sample)]
      }
      exprs5 <- exprs3
    } else{
      exprs5 <- data.frame(exprs3, Group = exprs3$Sample)
    }
  )

  # Aggregate data to group level by sample_aggregate_func.
  exprs6 <- data.table(exprs5)[, Exprs:=sample_aggregate_func(aggr_gene), by=list(Display, Group)] %>%
    setkey(Display, Group) %>%
    unique() %>%
    as.data.frame()

  ## Sort data.
  if(!is.null(sort_rows_by)){
    # create a new variable to split data on from sort_rows_by.
    avg_exprs <- tryCatch(expr = {
      wh <- strsplit(sort_rows_by,";") %>% unlist() %>%
      sapply(.,strsplit,split = "=") %>%
      sapply(.,function(x){
        paste0(x[1],"=='",x[2],"'")}) %>% paste(.,collapse = " & ")

      avg_exprs <- exprs6 %>%
        mutate(sort_by = eval(parse(text = paste0("ifelse(",wh,",'x','y')"))))
    },error = function(e){stop("Please ensure 'sort_rows_by' follows correct syntax and the variables are included in either 'group_by' or 'facet_by'.")})

    # Apply sort_rows_func for every unique Display value, grouped according to sort_rows_by.
    avg_exprs <- tryCatch(expr = {
      avg_exprs %>%
      select(Exprs, sort_by, Display) %>%
      group_by(Display,sort_by) %>%
      summarise(value = list(Exprs)) %>%
      spread(sort_by, value) %>%
      group_by(Display) %>%
      mutate(score = sort_rows_func(unlist(x),unlist(y))) %>%
      select(-x,-y) %>%
      arrange(desc(score))
    },error = function(e){stop("Please ensure 'sort_rows_func' takes two input vectors and outputs one value, used for ordering.")})
  } else {
    avg_exprs <- group_by(exprs6, Display) %>%
      summarise(score = sort_rows_func(Exprs)) %>%
      arrange(desc(score))
  }

  ## Subset to X genes, either by order or by a character vector.
  if (is.numeric(gene_show)){
    if (gene_show > nrow(avg_exprs)){
      gene_show <- nrow(avg_exprs)
    }
    exprs7 <- filter(exprs6, Display %in% unique(avg_exprs$Display)[1:gene_show]) %>%
      mutate(Display = {factor(Display,levels = rev(unique(avg_exprs$Display)))})
  }
  if (!is.numeric(gene_show)){
    if (all(gene_show != "all")){
      exprs7 <- filter(exprs6, Display %in% gene_show) %>%
        mutate(Display = {factor(Display,levels = rev(gene_show))})
    }
    if (all(gene_show == "all")){
      gene_show <- nrow(avg_exprs)
      exprs7 <- filter(exprs6, Display %in% avg_exprs$Display[1:gene_show])
    }
  }

  ## Define the output
  ## Make a heatmap style plot
    heatmap <- ggplot(exprs7, aes_string(x = "Group", y = "Display", label = formatC("Exprs", format = "f", digits = 1))) +
      geom_tile(aes(fill = Exprs), colour = "white", size = 0.5) +
      theme(axis.text.y = element_text(size = 12, color = "black", vjust = 0.4),
            axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, angle = 90, hjust = 1),
            axis.title = element_blank(),
            text = element_text(size = 8, color = "black"),
            axis.line = element_blank(),
            #axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(1,1,1,1), "mm"),
            title = element_text(size = 8),
            panel.background = element_blank())

    ## Get colorpalette for colorscale or set default
    if (!is.null(color_vector)){
      color.pal = color_vector
    } else {
      color.pal = rev(RColorBrewer::brewer.pal(3, "RdBu"))
    }

    if (plot_values == TRUE){
      exprs8 <- exprs7
      exprs8$Exprs <- round(exprs8$Exprs, round)
      heatmap <- heatmap + geom_text(data = exprs8, size = plot_values_size, colour = "grey10", check_overlap = TRUE) +
        theme(legend.position = "none")
    }
    # Colors of heatmap.
    heatmap <- heatmap + scale_fill_gradientn(
      colours  = color.pal,
      trans    = plot_colorscale,
      na.value = "#67A9CF",
      oob      = scales::squish)

    if(!is.null(facet_by)){
      if(length(ogroup) > 1){
        heatmap$data$Group <- apply(heatmap$data[,ogroup], 1, paste, collapse = " ")
      } else{
        heatmap$data$Group <- heatmap$data[,ogroup]
      }

      if(plot_values == TRUE){
        if(length(ogroup) > 1){
          heatmap$layers[[2]]$data$Group <- apply(heatmap$layers[[2]]$data[,ogroup], 1, paste, collapse = " ")
        } else{
          heatmap$layers[[2]]$data$Group <- heatmap$layers[[2]]$data[,ogroup]
        }
      }
      heatmap <- heatmap + facet_grid(reformulate(facet_by), scales = "free_x", space = "free")
      heatmap <- heatmap + theme(strip.text = element_text(size = 10))
    }
    return(heatmap)
}
