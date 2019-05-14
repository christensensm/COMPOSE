#' PCA and Pearson correlation plots from counts table
#'
#' Basic comparison of CRISPR screen counts: PCA and Pearson correlation plots
#'
#' @param countsTable input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data
#' @param identifier1 string identifying column name of metadata with which to adjust color in PCA plot
#' @param identifier2 string identifying column name of metadata with which to adjust alpha in PCA plot
#' @param identifier3 string identifying column name of metadata with which to adjust shape in PCA plot
#' @param save logical - do you want to save the violin plot to pdf
#' @return ggplot object of the violin plot
#' @seealso \code{\link{ggbiplot}} which this function uses to plot PCAs
#' @seealso \code{\link{ggplot}} which this function uses to plot Pearson correlations
#' @export
#' @examples
#' y <- matrix(rnorm(100*9, mean = 10, sd = 1),100,9)
#' y[,1:3] <- y[,1:3] + 5
#' metadata <- data.frame(sample = paste0('sample.',1:9), batch = c("A","A","A","B","B","B","C","C","C"))
#'
#' count.pca(y, metadata, identifier1 = 'sample')
#' count.pca(y, metadata, identifier1 = 'sample', batch = T, batch.id = 'batch')
#' ...
#' @importFrom stats prcomp cor as.dist hclust as.dendrogram
#' @importFrom grDevices hcl pdf dev.off
#' @importFrom ggbiplot ggbiplot
#' @importFrom ggplot2 geom_point scale_fill_manual scale_shape_manual scale_alpha_manual ggtitle theme theme_bw geom_tile aes guide_legend element_text element_blank
#' @importFrom reshape2 melt
#' @importFrom edgeR cpm
#' @importFrom ggdendro ggdendrogram
#' @import viridis
#' @importFrom gplots heatmap.2
#' @importFrom limma removeBatchEffect

# Create function to compute principal components and plot with pearson correlation
count.pca <- function(countsTable, metadata, identifier1 = NULL, identifier2 = NULL, identifier3 = NULL,
                      batch = FALSE, batch.id = NULL, save = F) {
    if (ncol(countsTable) != nrow(metadata))
        stop("Differing number of samples in count matrix and metadata table")
    if (!all(colnames(countsTable) == metadata[, 1]))
        stop("Please make sure your count matrix column names are the same as your metadata sample IDs (first column)")
    countsTable <- edgeR::cpm(countsTable)
    if(batch == TRUE){
      batch.col <- which(colnames(metadata) == batch.id)
      countsTable <- limma::removeBatchEffect(countsTable, batch = metadata[,batch.col])
    }
    df.pca = stats::prcomp(t(countsTable))  #create prcomp object
    # identify the column in metadata for each identifier
    identifier1.col = which(colnames(metadata) == identifier1)
    identifier2.col = which(colnames(metadata) == identifier2)
    identifier3.col = which(colnames(metadata) == identifier3)
    # set colors and shapes identifier1 = color identifier2 = alpha identifier3 = shape
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    if(is.null(identifier1))
      stop("Please select an identifier column name in the metadata for 'identifier1'")
    fillcol = gg_color_hue(length(unique(metadata[, identifier1])))
    names(fillcol) <- unique(metadata[, identifier1])
    shapes = 21:25
    # Determine number of identifiers (must be at least one identifier)
    if (!is.null(identifier2) & !is.null(identifier3)) {
        fill.id <- factor(metadata[, identifier1.col], levels = unique(metadata[, identifier1.col]))
        alpha.id <- factor(metadata[, identifier2.col], levels = unique(metadata[, identifier2.col]))
        shape.id <- factor(metadata[, identifier3.col])
        # Plot PC1 vs PC2
        p1 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 2), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC1 vs PC2") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p2 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC1 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC2 vs PC3
        p3 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(2, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC2 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Determine pearson correlations, reorder by hierarchical clustering, and plot
        datcor <- stats::cor(countsTable)
        # Plot all four graphs
        if (save) {
          grDevices::pdf("PCA_PC1vPC2_normalized_counts.pdf", width = 8, height = 8)
          print(p1)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC1vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p2)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC2vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p3)
          grDevices::dev.off()
          grDevices::pdf("Pearson_cor_normalized_counts.pdf", width = 8, height = 8)
          print(corheat)
          grDevices::dev.off()
        }
        print(p1)
        print(p2)
        print(p3)
        gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
    }
    if (!is.null(identifier2) & is.null(identifier3)) {
        fill.id <- factor(metadata[, identifier1.col], levels = unique(metadata[, identifier1.col]))
        alpha.id <- factor(metadata[, identifier2.col], levels = unique(metadata[, identifier2.col]))
        # Plot PC1 vs PC2
        fill.id <- factor(metadata[, identifier1.col], levels = unique(metadata[, identifier1.col]))
        alpha.id <- factor(metadata[, identifier2.col], levels = unique(metadata[, identifier2.col]))
        p1 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 2), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::ggtitle("PC1 vs PC2") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p2 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::ggtitle("PC1 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p3 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(2, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, alpha = alpha.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_alpha_manual(name = identifier2, values = seq(0, 1, 1/length(levels(alpha.id)))[-1], labels = levels(alpha.id)) +
          ggplot2::ggtitle("PC2 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Determine pearson correlations, reorder by hierarchical clustering, and plot
        datcor <- stats::cor(countsTable)
        # Plot all four graphs
        if (save) {
          grDevices::pdf("PCA_PC1vPC2_normalized_counts.pdf", width = 8, height = 8)
          print(p1)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC1vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p2)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC2vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p3)
          grDevices::dev.off()
          grDevices::pdf("Pearson_cor_normalized_counts.pdf", width = 8, height = 8)
          print(corheat)
          grDevices::dev.off()
        }
        print(p1)
        print(p2)
        print(p3)
        gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
    }
    if (is.null(identifier2) & !is.null(identifier3)) {
        fill.id <- factor(metadata[, identifier1.col], levels = unique(metadata[, identifier1.col]))
        shape.id <- factor(metadata[, identifier3.col])
        # Plot PC1 vs PC2
        p1 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 2), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC1 vs PC2") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p2 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC1 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p3 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(2, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, ggplot2::aes(fill = fill.id, shape = shape.id)) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = levels(fill.id), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::scale_shape_manual(name = identifier3, values = unique(shapes[as.numeric(shape.id)]), labels = levels(shape.id)) +
          ggplot2::ggtitle("PC2 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Determine pearson correlations, reorder by hierarchical clustering, and plot
        datcor <- stats::cor(countsTable)
        # Plot all four graphs
        if (save) {
          grDevices::pdf("PCA_PC1vPC2_normalized_counts.pdf", width = 8, height = 8)
          print(p1)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC1vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p2)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC2vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p3)
          grDevices::dev.off()
          grDevices::pdf("Pearson_cor_normalized_counts.pdf", width = 8, height = 8)
          print(corheat)
          grDevices::dev.off()
        }
        print(p1)
        print(p2)
        print(p3)
        gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
    }
    if (is.null(identifier2) & is.null(identifier3)) {
        # Plot PC1 vs PC2
        p1 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 2), alpha = 0) +
          ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = as.factor(metadata[, identifier1.col]))) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = unique(as.factor(metadata[, identifier1.col])), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::ggtitle("PC1 vs PC2") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC1 vs PC3
        p2 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(1, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = as.factor(metadata[, identifier1.col]))) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = unique(as.factor(metadata[, identifier1.col])), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::ggtitle("PC1 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Plot PC2 vs PC3
        p3 <- ggbiplot::ggbiplot(df.pca, pc.biplot = TRUE, var.axes = FALSE, choices = c(2, 3), alpha = 0) +
          ggplot2::geom_point(size = 4, shape = 21, ggplot2::aes(fill = as.factor(metadata[, identifier1.col]))) +
          ggplot2::scale_fill_manual(name = identifier1, values = fillcol, labels = unique(as.factor(metadata[, identifier1.col])), guide = ggplot2::guide_legend(override.aes = list(shape = 21))) +
          ggplot2::ggtitle("PC2 vs PC3") +
          ggplot2::theme_bw() +
          ggplot2::theme(aspect.ratio = 1)
        # Determine pearson correlations, reorder by hierarchical clustering, and plot
        datcor <- stats::cor(countsTable)
        # Plot all four graphs
        if (save) {
          grDevices::pdf("PCA_PC1vPC2_normalized_counts.pdf", width = 8, height = 8)
          print(p1)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC1vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p2)
          grDevices::dev.off()
          grDevices::pdf("PCA_PC2vPC3_normalized_counts.pdf", width = 8, height = 8)
          print(p3)
          grDevices::dev.off()
          grDevices::pdf("Pearson_cor_normalized_counts.pdf", width = 8, height = 8)
          print(corheat)
          grDevices::dev.off()
        }
        print(p1)
        print(p2)
        print(p3)
        gplots::heatmap.2(datcor, trace = 'none', col = viridis::viridis(100), margins=c(10,10), srtCol = 45)
    }
}

