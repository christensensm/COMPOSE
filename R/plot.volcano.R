#' Volcano plot for sgRNA fold-change results
#'
#' Visualization of fold-changes and p-values for sgRNAs
#'
#' @param results a calc.DESeq2.L2FC object
#' @param sig.cutoff adjusted p-value cutoff for sgRNA fold-changes
#' @param print logical - do you want to print the  plot to plots
#' @param save logical - do you want to save the  plot to pdf
#' @seealso \code{\link{ggplot}, \link{ggrepel}} which this function uses to plot
#' @export
#' @examples
#'
#' ...
#' @importFrom ggplot2 ggplot aes theme_bw geom_point geom_hline scale_size_continuous scale_color_gradient xlab ylab ggtitle
#' @importFrom ggrepel geom_text_repel


# Create function to plot volcano plot visualizing sgRNA fold-changes and p-values
plot_volcano <- function(results, print = T, save = F, sig.cutoff = 0.05, contrast = NULL, top = NULL) {
  if(is.null(contrast)){
    cat(paste(seq(1:ncol(results$L2FC)),colnames(results$L2FC), sep = '. '),sep="\n")
    contrast <- as.numeric(readline(prompt="Enter the column number for which comparison you'd like to plot: "))
  }
  if(is.null(top))
    top <- as.numeric(readline(prompt="Enter the number of top genes to label: "))
  volcano <- data.frame(l2fc = results$L2FC[,contrast], padj = results$padj[,contrast])
  volcano <- extract.genefromguide(volcano, neg.controls)
  volcano <- volcano[order(volcano$padj),]
  volcano$show.in.plot <- c(rep(1,times = top),rep(0,times = nrow(volcano)-top))
  volcano$show.in.plot[volcano$show.in.plot==1] <- volcano$gene.name[volcano$show.in.plot==1]
  volcano$show.in.plot[volcano$show.in.plot==0] <- ""
  p <- ggplot2::ggplot(volcano, ggplot2::aes(x = l2fc, y = -log10(padj), color = padj, size = -log10(padj))) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept=-log10(sig.cutoff), linetype="dashed", color = "red") +
    ggplot2::scale_size_continuous(guide = F) +
    ggplot2::scale_color_gradient(low = "#56B1F7", high = "#132B43") +
    ggplot2::theme_bw() +
    ggplot2::xlab("log2(Fold-Change)") + ggplot2::ylab("-log10(p-value)") +
    ggplot2::ggtitle(colnames(results$L2FC)[contrast]) +
    ggrepel::geom_text_repel(
      data = subset(volcano, show.in.plot!=""),
      size = 4,
      color = 'black',
      ggplot2::aes(x = subset(volcano, show.in.plot!="")$l2fc,
                   y = -log10(subset(volcano, show.in.plot!="")$padj),
                   label = subset(volcano, show.in.plot!="")$show.in.plot))
  if(print) print(p)
  if(save){
    grDevices::dev.copy(pdf,paste0('volcano_plot_sgRNAs_',colnames(results$L2FC)[contrast],'.pdf'),height = 8, width = 8)
    grDevices::dev.off()
  }
}
