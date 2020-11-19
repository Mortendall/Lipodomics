library(tidyverse)
library(lipidr)
library(ggplot2)
library(limma)
library(vroom)
library(openxlsx)
library(here)
library(fs)



#' Lipodomics data loader
#'
#' @param file_type the name of the count file
#'
#' @return a trimmed count matrix

count_file_loader <- function(file_type){
  count_file <- fs::dir_ls(here("data-raw/"),
                           regexp = file_type,
                           recurse = TRUE)
  count_data <- read.xlsx(count_file)
  rownames(count_data) <- count_data$`mol%`
  count_data <- count_data %>% dplyr::select(-`mol%`)
  return(count_data)
}

#' Load metadata
#'
#' @param file_name the name of the metadata file
#'
#' @return metadata file

load_metadata <- function(file_name) {
  data_file <- fs::dir_ls(here::here("data-raw/"),
                          regexp = file_name,
                          recurse = T)
  metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
  return(metadata)
}

#' Quality control generator
#'
#' @param count_matrix a count matrix generated through the count_matrix function
#' @param setup setup data.frame
#'
#' @return

Quality_control_plots <- function(count_matrix, setup) {
  group <- as.matrix(setup$Group)
  RNAseq <- edgeR::DGEList(counts = count_matrix, group = group)

  pD <-
    reshape2::melt(cpm(RNAseq, normalized.lib.sizes = TRUE, log = TRUE))
  p <- ggplot(pD, aes(value)) +
    geom_density() +
    facet_wrap( ~ Var2)
  dir.create(here::here("data/figures/QCplots"), showWarnings = F)
  ggplot2::ggsave(
    p,
    filename = here("data/figures/QCplots/Density_plot.pdf"),
    width = 8,
    height = 8,
    units = "cm",
    scale = 1
  )

  #Create mdPlots

  oldpar <- par()$mfrow
  pdf(file.path(here("data/figures/QCplots"), "beforeFiltering_MD.pdf"), width = 4, height = 4)
  par(mfrow = c(2, 2))
  for (i in seq_len(ncol(RNAseq))) {
    plotMD(RNAseq, column = i)
    abline(h = 0)
  }
  par(mfrow = oldpar)
  dev.off()

  #create Possion heatmap
  poisd <- PoissonDistance(t(RNAseq$counts))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(cpm(RNAseq))
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  heatmap <- pheatmap(samplePoisDistMatrix,  clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd, col=colors)
  ggsave(heatmap, filename = here("data/figures/QCplots/Poisson_heatmap.png"),
         width = 12,
         height = 12,
         units = "cm",
         scale = 2.5)
  #crete mdsPlots
  mdsData <- plotMDS(RNAseq, ndim = 3, plot = FALSE)
  mdsData <-
    mdsData$cmdscale.out %>% data.table(keep.rownames = TRUE) %>%
    mutate(ID = rownames(RNAseq$samples)) %>%
    dplyr::select(-rn) %>%
    mutate(Group = setup$Group)

  setnames(mdsData,
           c("V1", "V2", "V3", "ID", "Group"),
           c("dim1", "dim2", "dim3", "ID", "Group"))
  plotMDS(RNAseq, ndim = 3)
  pBase <-
    ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()
  pBase2 <-
    ggplot(mdsData, aes(x = dim1, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pBase3 <-
    ggplot(mdsData, aes(x = dim2, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pdf(file.path(here("data/figures/QCplots"), "MDSplots.pdf"), width = 4, height = 4)
  par(mfrow = c(1, 1))
  plot(pBase)
  plot(pBase2)
  plot(pBase3)
  plotMDS(RNAseq, ndim = 3)
  par(mfrow = oldpar)
  dev.off()



  print("All your plots can be found in the Figures/QCplots folder")
}

#' Differentially expressed gene analysis
#'
#' @param data a count matrix
#' @param meta setup data
#' @param proteome_key
#'
#' @return a DEG analysis
meta <- setup
data <- count_matrix
all(setup$ID==colnames(count_matrix))
DEG_analysis <- function(data, meta){
  design <- stats::model.matrix( ~0+Group, meta)
  colnames(design) <- str_remove_all(colnames(design), "Group")
  fit <- limma::lmFit(data, design = design, method = "robust")
  fit<- eBayes(fit)

  cont.matrix <- makeContrasts(
    Treatment_in_WT = WT_NR - WT_Control,
    Treatment_in_KO = KO_NR - KO_Control,
    KO_in_Control = KO_Control - WT_Control,
    KO_in_NR = KO_NR - WT_NR,
    Interaction = (KO_NR - KO_Control) - (WT_NR - WT_Control),
    levels = design
  )

  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  resultTables <- list(
    Treatment_in_WT = topTable(fit2, coef = "Treatment_in_WT", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Treatment_in_KO = topTable(fit2, coef = "Treatment_in_KO", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_Control = topTable(fit2, coef = "KO_in_Control", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_NR = topTable(fit2, coef = "KO_in_NR", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
  )
  View(resultTables)
  write.xlsx(resultTables, file = here("data/limma_results.xlsx"))
  return(NASH_vs_con)
}


#' Limma result data loader
#'
#' @param file_name name of the limma_results file
#'
#' @return limma analysis

load_limma_data <- function(file_name) {
  data_file <- fs::dir_ls(here("data/"),
                          regexp = file_name,
                          recurse = TRUE)
  limma_analysis <-read.xlsx(data_file)
  return(limma_analysis)
}

#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list results table generated from limma. Must be data.table and contain a SYMBOL annotation column
#'
#' @return a list containing enrichresults for each element in the results file list


goAnalysis <- function(result_list){
  result_list <- as.data.table(result_list)
  bg_list <- clusterProfiler::bitr(
    result_list$Gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db",
    drop = T
  )

  sig_list<- result_list %>%
    dplyr::filter(adj.P.Val<0.05)

  eg <- clusterProfiler::bitr(
    sig_list$Gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db",
    drop = T
  )
  goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                         universe = bg_list$ENTREZID,
                                         OrgDb = org.Hs.eg.db,
                                         ont = "BP")

  return(goResults)
}
