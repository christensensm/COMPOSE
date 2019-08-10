#' Upset plot for gene score results
#'
#' Visualization of gene scores and overlaps
#'
#' @param counts input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data including batch, treatments, etc.
#' @param sample.id string identifying column name in metadata containing sample IDs
#' @param print logical - do you want to print the  plot to plots
#' @param save logical - save the  plot to pdf
#' @param identifier1 string identifying column name of metadata with which to adjust color in PCA plot
#' @param identifier2 string identifying column name of metadata with which to adjust size in PCA plot
#' @param identifier3 string identifying column name of metadata with which to adjust shape in PCA plot
#' @param batch logical - correct for batch effects (requires a batch.id input)
#' @param batch.id numerical - column of metadata that identifies the batch effect to remove
#' @param cpm logical - do you want to normalize the counts by sequencing depth
#' @param transform logical - do you want to log10 transform the counts
#' @param thresh numerical cutoff for low counts
#' @param minsample.ids minimum number of samples for which counts must be present per gene
#' @param meta.idnum vector containing the column numbers in metadata that represent (1) Cell line, (2) replicates, and (3) condition
#' @param cutoff gene score cutoff for comparisons between contrasts (default = 1)
#' @param sig.cutoff numeric - specified p-value cut-off
#' @param verbose TRUE/FALSE (default = TRUE)
#' @param top.labels numerical - number of top guides to label (by p-value)
#' @param neg.controls input vector containing gRNA ids for negative controls
#' @param essential.genes input vector containing gRNA ids for essential genes
#' @param split defines the regex used to split the guide RNA id to extract the gene name
#' @param pathway.analysis (TRUE/FALSE) continue with pathway analysis for all contrasts
#' @param gene.db genome database from which to convert gene symbols to other ids
#' @param list.return logical - return a list containing gene names for every intersection
#' @seealso \code{\link{COMPOSE}} the COMPOSE package
#' @export
#' @examples
#'
#' ...
#' @importFrom grDevices dev.copy dev.off pdf


# Create function to plot volcano plot visualizing sgRNA fold-changes and p-values
conduct <- function(counts = NULL, metadata = NULL, sample.id = "SampleID", print = T, save = F,
                    identifier1 = NULL, identifier2 = NULL, identifier3 = NULL, batch = F, batch.id = NULL,
                    transform = T, thresh = 1, minsample.ids = 2,
                    meta.idnum = NULL,
                    cutoff = 0.5, sig.cutoff = 0.05, verbose = T,
                    top.labels = 30, neg.controls = NULL, essential.genes = NULL, split = "_",
                    list.return = T, pathway.analysis = F, gene.db = NULL) {
  grDevices::pdf(paste0("orchestra_",Sys.Date(),".pdf"))
  #Plot sequencing depth per sample
  COMPOSE::calc.seq.depth(counts, metadata, sample.id = sample.id, print = print, save = save)

  #Violin plot visualizing read distribution
  COMPOSE::count.violin(counts, metadata, identifier1 = identifier1,
               identifier2=identifier2, transform = transform, save = save)

  #Compute principal components and plot with pearson correlation
  COMPOSE::count.pca(counts, metadata, identifier1=identifier1,
            identifier2=identifier2, identifier3=identifier3, batch = batch,
            batch.id = batch.id, save = save)

  counts.filtered <- COMPOSE::filterCounts(counts, thresh = thresh, minsample.ids = minsample.ids)

  results <- COMPOSE::calc.DESeq2.L2FC(counts.filtered, metadata = metadata, meta.idnum = meta.idnum,
                              include.batch = batch,
                              p.cutoff = sig.cutoff, save = save, verbose = verbose)

  COMPOSE::plot_sigguides(results, print = print, cutoff = cutoff, sig.cutoff = sig.cutoff, save = save)
  for(i in 1:ncol(results)){
    COMPOSE::plot_volcano(results, print = print, save = save, contrast = i, top = top.labels)
  }

  genescore <- COMPOSE::guides2genes(results$L2FC, neg.controls = neg.controls, essential.genes = essential.genes,
                            save.plots = save, save.tables = save, print.plots = print, split = split)

  COMPOSE::L2FC.violin(genescore, input = 'Gene Scores', variable = colnames(metadata)[meta.idnum[3]],
              save.plot = save, print.plot = print)

  intersections <- COMPOSE::plot_upset(genescore, list.return = list.return, save = save, cutoff = cutoff)

  dev.off()
  if(pathway.analysis){
    library(gene.db)
    pdf('pathway.analyses',Sys.Date(),'.pdf')
    COMPOSE::reactomepa.enrich(genescore[!rownames(genescore) %in% essential.genes,],
                               genome.db = gene.db, showCategory = top.labels,
                               cutoff = cutoff, pvalue = sig.cutoff, save = save)
    for(i in 1:ncol(genescore)){
      COMPOSE::reactomepa.enrich(genescore[!rownames(genescore) %in% essential.genes,i,drop=FALSE],
                                 genome.db = gene.db, showCategory = top.labels,
                                 cutoff = cutoff, pvalue = sig.cutoff, save = save)
    }
    dev.off()
  }
  results$genescore <- genescore
  return(results)
}


