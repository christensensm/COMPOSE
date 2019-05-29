#' Upset plot for gene score results
#'
#' Visualization of gene scores and overlaps
#'
#' @param genescore a guides2genes object
#' @param cutoff gene score cutoff for comparisons between contrasts (default = 1)
#' @param save logical - save the  plot to pdf
#' @param list.return logical - return a list containing gene names for every intersection
#' @seealso \code{\link{UpSetR}, \link{Vennerable}} which this function uses to plot
#' @export
#' @examples
#'
#' ...
#' @importFrom grDevices dev.copy dev.off
#' @importFrom UpSetR upset
#' @importFrom Vennerable Venn


# Create function to plot volcano plot visualizing sgRNA fold-changes and p-values
plot_upset <- function(genescore, save = FALSE, cutoff = 1, list.return = FALSE) {
  #Positively enriched genes (score >= cutoff)
  pos.enrich.upset = data.frame(Gene = rownames(genescore), genescore)
  pos.enrich.upset[,-1] = as.numeric(as.logical(pos.enrich.upset[,-1]>=cutoff))
  UpSetR::upset(pos.enrich.upset, order.by = 'freq', nsets = ncol(pos.enrich.upset),
                empty.intersections = 'on')
  if(save){
    grDevices::dev.copy(pdf,paste0('upset_plot_enriched_genes.pdf'),height = 8, width = 8)
    grDevices::dev.off()
  }
  if(list.return){
    genes.up <- list()
    cutoff=1
    for(i in 1:ncol(genescore)){
      genes.up[[i]] <- rownames(genescore[genescore[,i] >= cutoff,])
    }
    names(genes.up) <- colnames(genescore)
    temp <- Vennerable::Venn(genes.up)
    temp@IntersectionSets <- temp@IntersectionSets[order(temp@IndicatorWeight[,ncol(temp@IndicatorWeight)],decreasing = T)]
    temp@IndicatorWeight <- temp@IndicatorWeight[order(temp@IndicatorWeight[,ncol(temp@IndicatorWeight)],decreasing = T),]
    temp.names <- names(temp@IntersectionSets)
    for(i in 1:length(temp.names)){
      if(sum(as.numeric(temp.names[i]))==0)
        temp.names[i] <- 'none'
      else
        temp.names[i] <- paste(colnames(temp@IndicatorWeight)[-ncol(temp@IndicatorWeight)][temp@IndicatorWeight[i,-ncol(temp@IndicatorWeight)]==1],
                               collapse = ".")
    }
    names(temp@IntersectionSets) <- paste0(temp.names,"_enriched")
    temp@IntersectionSets <- temp@IntersectionSets[names(temp@IntersectionSets)!='none_enriched']
  }

  #Negative enriched genes (score <= -1)
  neg.enrich.upset = data.frame(Gene = rownames(genescore), genescore)
  neg.enrich.upset[,-1] = as.numeric(as.logical(neg.enrich.upset[,-1]<=-cutoff))
  UpSetR::upset(neg.enrich.upset, order.by = 'freq', nsets = ncol(neg.enrich.upset),
                empty.intersections = 'on')
  if(save){
    grDevices::dev.copy(pdf,paste0('upset_plot_depleted_genes.pdf'),height = 8, width = 8)
    grDevices::dev.off()
  }
  if(list.return){
    genes.down <- list()
    cutoff=1
    for(i in 1:ncol(genescore)){
      genes.down[[i]] <- rownames(genescore[genescore[,i] <= -cutoff,])
    }
    names(genes.down) <- colnames(genescore)
    temp2 <- Vennerable::Venn(genes.down)
    temp2@IntersectionSets <- temp2@IntersectionSets[order(temp2@IndicatorWeight[,ncol(temp2@IndicatorWeight)],decreasing = T)]
    temp2@IndicatorWeight <- temp2@IndicatorWeight[order(temp2@IndicatorWeight[,ncol(temp2@IndicatorWeight)],decreasing = T),]
    temp2.names <- names(temp2@IntersectionSets)
    for(i in 1:length(temp2.names)){
      if(sum(as.numeric(temp2.names[i]))==0)
        temp2.names[i] <- 'none'
      else
        temp2.names[i] <- paste(colnames(temp2@IndicatorWeight)[-ncol(temp2@IndicatorWeight)][temp2@IndicatorWeight[i,-ncol(temp2@IndicatorWeight)]==1],
                                collapse = ".")
    }
    names(temp2@IntersectionSets) <- paste0(temp2.names,"_depleted")
    temp2@IntersectionSets <- temp2@IntersectionSets[names(temp2@IntersectionSets)!='none_depleted']

    gene.intersection.list <- c(temp@IntersectionSets, temp2@IntersectionSets)
    return(gene.intersection.list)
  }
}


