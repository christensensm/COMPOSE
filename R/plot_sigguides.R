#' Bidirectional barplot for sgRNA fold-change results
#'
#' Visualize the number of significantly enriched/depleted sgRNAs
#'
#' @param results a calc.DESeq2.L2FC object
#' @param sig.cutoff adjusted p-value cutoff for sgRNA fold-changes
#' @param fc.cutoff log2 fold-change cutoff for sgRNAs
#' @param print logical - do you want to print the  plot to plots
#' @param save logical - do you want to save the  plot to pdf
#' @seealso \code{\link{ggplot}} which this function uses to plot
#' @export
#' @examples
#'
#' ...
#' @importFrom ggplot2 ggplot aes theme_bw geom_point geom_hline scale_size_continuous scale_color_gradient xlab ylab ggtitle
#' @importFrom ggrepel geom_text_repel


# Create function to plot volcano plot visualizing sgRNA fold-changes and p-values
plot_sigguides <- function(results, print = T, save = F, sig.cutoff = 0.05, cutoff = 1) {
  numsig.guides <- data.frame(names = rep(colnames(results$L2FC)[1], 2),
                              Direction = factor(c("Depleted","Enriched"),levels=c("Enriched","Depleted")),
                              sig.guides = c(sum(results$padj[,1]<=sig.cutoff & results$L2FC[,1] <= -cutoff),
                                             sum(results$padj[,1]<=sig.cutoff & results$L2FC[,1] >= cutoff)))
  for(i in 2:ncol(results$L2FC)){
    temp.dat <- data.frame(names = rep(colnames(results$L2FC)[i], 2),
                           Direction = c("Depleted","Enriched"),
                           sig.guides = c(sum(results$padj[,i]<=sig.cutoff & results$L2FC[,i] <= -cutoff),
                                          sum(results$padj[,i]<=sig.cutoff & results$L2FC[,i] >= cutoff)))
    numsig.guides <- rbind(numsig.guides, temp.dat)
  }
  numsig.guides$sig.guides <- ifelse(numsig.guides$Direction == "Depleted",
                                     -1 * numsig.guides$sig.guides,
                                     numsig.guides$sig.guides)
  cols = c("Depleted" = "orangered","Enriched" = "forestgreen")
  p <- ggplot2::ggplot(numsig.guides, ggplot2::aes(fill=Direction,
                                              y=sig.guides, x=names)) +
    ggplot2::geom_bar(stat="identity") + theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = margin(2,2,10,80), axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::ylab("Number of significant guides")


  if(print) print(p)
  if(save){
    grDevices::dev.copy(pdf,'significant_sgRNAs_barplot.pdf',height = 8, width = 8)
    grDevices::dev.off()
  }
}
