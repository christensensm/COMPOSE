#' Count Distribution
#'
#' Violin plot visualizing read distribution
#'
#' @param countsTable input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data
#' @param identifier1 string identifying column name of metadata with which to adjust color in violin plot
#' @param identifier2 string identifying column name of metadata with which to adjust alpha in violin plot
#' @param cpm logical - do you want to normalize the counts by sequencing depth
#' @param transform logical - do you want to log10 transform the counts
#' @param save logical - do you want to save the violin plot to pdf
#' @return ggplot object of the violin plot
#' @seealso \code{\link{ggplot}} which this function uses to plot
#' @export
#' @examples
#' y <- matrix(rnorm(100*9, mean = 10, sd = 1),100,9)
#' colnames(y) <- paste0('sample.',1:9)
#' metadata <- data.frame(sample = paste0('sample.',1:9), batch = c("A","A","A","B","B","B","C","C","C"),
#' time = rep(c(1,2,3),times = 3))
#' count.violin(y, metadata, identifier1 = 'sample', identifier2='time', transform = T)
#' ...
#' @importFrom ggplot2 ggplot aes theme_bw geom_violin geom_boxplot theme element_text element_blank margin
#' @importFrom grDevices dev.copy dev.off
#' @importFrom edgeR cpm


# Create function to plot violin plot visualizing read distribution
count.violin <- function(countsTable, metadata, identifier1 = "sample", identifier2 = NULL, transform = FALSE, cpm = FALSE, save = F) {
    if (ncol(countsTable) != nrow(metadata))
        stop("Differing number of samples in count matrix and metadata table")
    if (!all(colnames(countsTable) == metadata[, 1]))
        stop("Please make sure your count matrix column names are the same as your metadata sample IDs (first column)")
    if(cpm)
      countsTable <- edgeR::cpm(countsTable)
    # create sample id vector
    sample.ids <- colnames(countsTable)
    for (i in 1:length(sample.ids)) {
        if (i == 1)
            sample.id.dat <- c(rep(sample.ids[i], times = nrow(countsTable)))
        else
          sample.id.dat <- c(sample.id.dat, rep(sample.ids[i], times = nrow(countsTable)))
    }
    # create counts vector
    for (i in 1:ncol(countsTable)) {
        if (i == 1)
            count.dat <- c(countsTable[, i])
        else
          count.dat <- c(count.dat, countsTable[, i])
    }
    # create identifier 1 vector
    identifier1 = which(colnames(metadata) == identifier1)
    if (length(identifier1) == 0)
        stop("Identifier 1 does not exist in metadata")
    for (i in 1:length(metadata[, identifier1])) {
        if (i == 1)
            identifier1.dat <- c(rep(as.character(metadata[, identifier1][i]), times = nrow(countsTable))) else identifier1.dat <- c(identifier1.dat, rep(as.character(metadata[, identifier1][i]), times = nrow(countsTable)))
    }
    # Run if identifier2 is specified
    if (!is.null(identifier2)) {
        # create identifier 2 vector
        identifier2 = which(colnames(metadata) == identifier2)
        if (length(identifier2) == 0)
            stop("Identifier 2 does not exist in metadata")
        for (i in 1:length(metadata[, identifier2])) {
            if (i == 1)
                identifier2.dat <- c(rep(as.character(metadata[, identifier2][i]), times = nrow(countsTable))) else identifier2.dat <- c(identifier2.dat, rep(as.character(metadata[, identifier2][i]), times = nrow(countsTable)))
        }
        # Create dataframe that contains columns for Log2FC and sample (plus identifier1 and 2)
        norm.dat <- data.frame(counts = count.dat,
                               sample.id = factor(sample.id.dat, levels = sample.ids),
                               identifier1 = identifier1.dat, identifier2 = identifier2.dat)
        # Create violin plot using ggplot
        if (transform == T)
            norm.dat$counts <- log10(norm.dat$counts + 1)
        # creates ggplot transforms to violin plot with colors based on identifier1, alpha based on identifier2 Boxplot added
        # transforms to violin plot with colors based on identifier1
        p <- ggplot2::ggplot(norm.dat, ggplot2::aes(x = factor(sample.id, levels = unique(norm.dat$sample.id)), y = counts)) +
          ggplot2::theme_bw() +
          ggplot2::geom_violin(ggplot2::aes(fill = identifier1,
            alpha = 1 - (as.numeric(identifier2)/length(unique(identifier2))))) +
          ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none", axis.title.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(2,2,10,80)) +
        print(p)  #plots
        if (save) {
          grDevices::dev.copy(pdf, "Read_Distribution_per_sample.pdf", height = 8, width = nrow(metadata))
          grDevices::dev.off()
        }
    }
    # Run if no identifier2 is specified
    if (is.null(identifier2)) {
        # Create dataframe that contains columns for Log2FC and sample (plus identifier1)
        norm.dat <- data.frame(counts = count.dat, sample.id = factor(sample.id.dat, levels = sample.ids),
            identifier1 = identifier1.dat)
        # Create violin plot using ggplot
        if (transform == T)
            norm.dat$counts <- log10(norm.dat$counts + 1)
        # creates ggplot transforms to violin plot with colors based on identifier1 Boxplot added
        p <- ggplot2::ggplot(norm.dat, aes(x = factor(sample.id, levels = unique(norm.dat$sample.id)),
                                           y = counts)) +
          ggplot2::theme_bw() +
          ggplot2::geom_violin(aes(fill = identifier1)) +
          ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
          ggplot2::ggtitle("Read Distribution per Sample") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         legend.position = "none", axis.title.x = ggplot2::element_blank(),
                         plot.margin = ggplot2::margin(2, 2, 10, 80)) +
        print(p)  #plots
        if (save) {
          grDevices::dev.copy(pdf, "Read_Distribution_per_sample.pdf", height = 8, width = nrow(metadata))
          grDevices::dev.off()
        }
    }
}
