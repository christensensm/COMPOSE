#' Log2 Fold-Change Distribution
#'
#' Violin plot visualizing log2 fold-change values
#'
#' @param L2FC input matrix containing log2 fold-change values
#' @param input string identifying what type of input L2FC is (this function can also be used for visualizing CRISPR scores and CRISPR screen scores)
#' @param variable string identifying what variable was used to calculate L2FC values (for plot title purposes only)
#' @param print.plot logical - do you want to print the violin plot
#' @param save.plot logical - do you want to save the violin plot to pdf
#' @return ggplot object of the violin plot
#' @seealso \code{\link{ggplot}} which this function uses to plot
#' @export
#' @examples
#' L2FC.violin(L2FC.table)

# Create function to plot violin plot visualizing L2FC (or other scores)
L2FC.violin <- function(L2FC, input = "L2FC", variable = "time", save.plot = F, print.plot = T) {
    # create sample vector
    for (i in 1:ncol(L2FC)) {
        if (i == 1)
            sample.dat <- c(rep(colnames(L2FC)[i], times = nrow(L2FC))) else sample.dat <- c(sample.dat, rep(colnames(L2FC)[i], times = nrow(L2FC)))
    }
    # create counts vector
    for (i in 1:ncol(L2FC)) {
        if (i == 1)
            count.dat <- c(L2FC[, i]) else count.dat <- c(count.dat, L2FC[, i])
    }
    # create dataframe for violin plot
    norm.dat <- data.frame(gRNAname = rownames(L2FC), counts = count.dat, sample = factor(sample.dat, levels = colnames(L2FC)))
    # Create violin plot using ggplot creates ggplot transforms to violin plot with colors based on samples Add boxplot
    # adjust labels, remove legend/axis titles
    p <- ggplot2::ggplot(norm.dat, ggplot2::aes(x = sample, y = counts)) +
      ggplot2::theme_bw() +
      ggplot2::geom_violin(ggplot2::aes(fill = sample)) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none", axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(2,2,10,80)) +
      ggplot2::ggtitle(paste0(input, " - ", variable, "  comparison"))
    if (save.plot) {
        grDevices::pdf(paste0(input, "_", variable, "_comparison", format(Sys.time(), "%Y%b%d"), ".pdf"), width = max(ncol(L2FC), 8),
            height = 8)
        print(p)  #plots
        grDevices::dev.off()
    }
    if (print.plot)
        p
}

