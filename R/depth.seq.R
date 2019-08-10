#' Sequencing depth
#'
#' bar plot visualizing sequencing depth
#'
#' @param countsTable input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data
#' @param sample.id string identifying column name in metadata containing sample IDs
#' @param print logical - do you want to print the  plot to plots
#' @param save logical - do you want to save the  plot to pdf
#' @seealso \code{\link{ggplot}} which this function uses to plot
#' @export
#' @examples
#' y <- matrix(rnorm(100*9, mean = 10, sd = 1),100,9)
#' colnames(y) <- paste0('sample.',1:9)
#' metadata <- data.frame(sample = paste0('sample.',1:9), batch = c("A","A","A","B","B","B","C","C","C"))
#' calc.seq.depth(y)
#' ...
#' @importFrom ggplot2 ggplot aes theme_bw geom_bar theme ggtitle element_text element_blank


# Create function to plot violin plot visualizing read distribution
calc.seq.depth <- function(countsTable, metadata, sample.id = 'SampleID', print = T, save = F) {
  if (ncol(countsTable) != nrow(metadata))
    stop("Differing number of samples in count matrix and metadata table")
  if (!all(colnames(countsTable) == metadata[, 1]))
    stop("Please make sure your count matrix column names are the same as your metadata sample IDs (first column)")
  samp.column <- which(colnames(metadata) == sample.id)
  seq.depth <- data.frame(sample = metadata[,samp.column], counts = colSums(countsTable))
  p <- ggplot2::ggplot(seq.depth, ggplot2::aes(x = factor(sample, levels = unique(seq.depth$sample)), y = counts)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = sample)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),
                   legend.position = 'none', axis.title.x = ggplot2::element_blank()) +
    ggplot2::coord_fixed(ratio = 20/(2*ncol(countsTable))) +
    ggplot2::ggtitle("Sequencing depth per sample")
  if(print) print(p)
  if(save){
    grDevices::dev.copy(pdf,'Sequencing_depth_per_sample.pdf',height = 8, width = nrow(metadata))
    grDevices::dev.off()
  }
}
