#' CRISPR Screen Scoring
#'
#' Creates array containing log2 fold-change input, CRISPR scores, and CRISPR screen scores (similar to z-scores)
#'
#' @param L2FC input matrix containing log2 fold-change data with gRNA ids as row names
#' @param neg.controls input vector containing gRNA ids for negative controls
#' @param essential.genes input vector containing gRNA ids for essential genes
#' @param save.plots save the plots to pdf (logical, default = FALSE)
#' @param print.plots print the plots (logical, default = TRUE)
#' @param save.tables save the CS and CSS tables (logical, default = FALSE)
#' @param split defines the regex used to split the guide RNA id to extract the gene name
#' @return array containing additional object of a gene score table
#' @seealso \code{\link{L2FC.violin}, \link{ggplot2}, \link{ggbiplot}}
#' @export
#' @examples
#' test.array <- create.CSS(L2FC.filtered, neg.controls = neg.controls,
#' essential.genes = essential.genes)
#' ...
#' @importFrom stats aggregate prcomp cor
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off hcl
#' @importFrom ggbiplot ggbiplot
#' @importFrom ggplot2 geom_point aes scale_fill_manual guide_legend ggtitle theme_bw theme geom_tile ggplot element_text element_blank


# Function to: 1) plot L2FC for negative controls and essential genes 2) calculate CRISPR score and screen score and plot
# for essential genes 3) plot PCA and correlation for CSS values 4) return array with L2FC, CS, CSS and FDR tables
guides2genes <- function(L2FC, neg.controls, essential.genes, save.plots = FALSE,
                       print.plots = TRUE, save.tables = FALSE, split = "_") {
  #Identify genes for each guide
  L2FC <- extract.genefromguide(L2FC, neg.controls = neg.controls, split = split)
  #Create gene score table
  genescore <- data.frame(gene = unique(L2FC$gene.name))
  for(i in 1:(ncol(L2FC)-1)){
    dat.tab <- data.frame(logFC = L2FC[,i], gene.name = L2FC[,ncol(L2FC)])
    fc.agg <- stats::aggregate(dat.tab$logFC, list(dat.tab$gene.name), mean)
    genescore[,i+1] <- fc.agg$x
    colnames(genescore)[i+1] <- colnames(L2FC)[i]
  }
  rownames(genescore) <- genescore[,1]
  genescore <- genescore[,-1]
  if(sum(is.na(genescore))>0){
    genescore <- na.omit(genescore)
    cat("Eliminating rows containing NAs", "\n")
  }
  # Save table
  if (save.tables) {
    utils::write.csv(genescore, file = paste0('CRISPR_GENE_SCORE', format(Sys.time(), "%Y%b%d"), ".csv"))
  }
  # Print and/or save plots of negative controls, essential gene dropout and CSS correction, and PCA/correlation
  if (save.plots | print.plots) {
    cat("Plotting", "\n")
    # Visualize L2FC of negative control gRNAs
    p1 <- L2FC.violin(L2FC[L2FC$gene.name %in% neg.controls, -ncol(L2FC)],
                      print.plot = T, save.plot = F) +
      ggtitle("L2FC - Negative Controls")
    # Visualize L2FC of essential gRNAs
    L2FC.essential <- L2FC[L2FC$gene.name %in% essential.genes, ]
    p2 <- L2FC.violin(L2FC.essential[, -ncol(L2FC.essential)],
                      print.plot = T, save.plot = F) +
      ggtitle("L2FC - 1210 Essential Genes")
    # Display CSS values for essential genes per sample as violin plots
    score.essential <- genescore[rownames(genescore) %in% essential.genes, ]
    p4 <- L2FC.violin(score.essential, print.plot = T, save.plot = F) +
      ggtitle("CRISPR gene score - 1210 Essential Genes")
    # Save plots to individual pdfs
    if (save.plots) {
      grDevices::pdf("Negative_controls.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p1)
      grDevices::dev.off()
      grDevices::pdf("L2FC_1210_essential_genes.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p2)
      grDevices::dev.off()
      grDevices::pdf("CRISPR_screen_score_1210_essential_genes.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p4)
      grDevices::dev.off()
    }
    # Plot all four graphs
    if (print.plots){
      print(p1)
      print(p2)
      print(p4)
    }
    # Calculate correlation and display correlation of CSS between all sample.ids, PCA as well
    df.pca = stats::prcomp(t(genescore))  #create prcomp object
    # set colors
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    fillcol = gg_color_hue(ncol(genescore))
    names(fillcol) <- colnames(genescore)
    # Plot PC1 vs PC2
    p1 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 2), alpha = 0) +
      ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = colnames(genescore))) +
      ggplot2::scale_fill_manual(name = "Comparisons", values = fillcol, labels = colnames(genescore), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
      ggplot2::ggtitle("PC1 vs PC2") +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
    # Plot PC1 vs PC3
    p2 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 3), alpha = 0) +
      ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = colnames(genescore))) +
      ggplot2::scale_fill_manual(name = "Comparisons", values = fillcol, labels = colnames(genescore), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
      ggplot2::ggtitle("PC1 vs PC3") +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
    # Plot PC2 vs PC3
    p3 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(2, 3), alpha = 0) +
      ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = colnames(genescore))) +
      ggplot2::scale_fill_manual(name = "Comparisons", values = fillcol, labels = colnames(genescore), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
      ggplot2::ggtitle("PC2 vs PC3") +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
    # Determine pearson correlations, reorder by hierarchical clustering, and plot
    datcor <- stats::cor(genescore)
    # Save and print plots
    if (save.plots) {
      grDevices::pdf("CRISPR_screen_score_PC1vPC2.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p1)
      grDevices::dev.off()
      grDevices::pdf("CRISPR_screen_score_PC1vPC3.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p2)
      grDevices::dev.off()
      grDevices::pdf("CRISPR_screen_score_PC2vPC3.pdf", height = 8, width = max(8, ncol(L2FC)))
      print(p3)
      grDevices::dev.off()
      grDevices::pdf("CRISPR_screen_score_pearson_correlation.pdf", height = 8, width = max(8, ncol(L2FC)))
      gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
      grDevices::dev.off()
    }
    print(p1)
    print(p2)
    print(p3)
    gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
  }
  cat("Done!")
  return(genescore)
}
