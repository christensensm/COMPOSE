#' Conduct: Single-function processing of CRISPR screen data
#'
#' This function ties together all the basics of the COMPOSE package to give you a one-stop shop for processing data.
#'
#' The required arguments to run the simplest analysis include counts, metadata, meta.idnum,
#' essential.genes, and neg.controls (if pathway.analysis = F).
#'
#' @param counts input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data including batch, treatments, etc.
#' @param sample.id string identifying column name in metadata containing sample IDs
#' @param save.plot logical - save the plot to pdf or print to screen (default saves everything to one pdf)
#' @param save.table logical - save the tables (default saves everything)
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
conduct <- function(counts = NULL, metadata = NULL, sample.id = "SampleID", save.plot = T, save.table = T,
                    identifier1 = NULL, identifier2 = NULL, identifier3 = NULL, batch = F, batch.id = NULL,
                    transform = T, cpm = T, thresh = 1, minsample.ids = 2,
                    meta.idnum = NULL,
                    cutoff = 0.5, sig.cutoff = 0.05, verbose = T,
                    top.labels = 30, neg.controls = NULL, essential.genes = NULL, split = "_",
                    list.return = T, pathway.analysis = F, gene.db = NULL) {
  if(save.plot)
    grDevices::pdf(paste0("COMPOSE_conduct_",Sys.Date(),".pdf"), width = max(8,ncol(counts)/2))
  #Plot sequencing depth per sample
  if(verbose)
    print("Plotting sequencing depth, read distribution, PCAs, and heatmaps")
  COMPOSE::calc.seq.depth(counts, metadata, sample.id = sample.id)

  #Violin plot visualizing read distribution
  COMPOSE::count.violin(counts, metadata, identifier1 = identifier1,
               identifier2=identifier2, transform = transform, cpm = cpm)

  #Compute principal components and plot with pearson correlation
  COMPOSE::count.pca(counts, metadata, identifier1=identifier1,
            identifier2=identifier2, identifier3=identifier3, batch = batch,
            batch.id = batch.id)
  if(verbose)
    print("Filtering counts")
  counts.filtered <- COMPOSE::filterCounts(counts, thresh = thresh, minsample.ids = minsample.ids)

  results <- COMPOSE::calc.DESeq2.L2FC(counts.filtered, metadata = metadata, meta.idnum = meta.idnum,
                              include.batch = batch,
                              p.cutoff = sig.cutoff, save = save.table, verbose = verbose)
  if(verbose)
    print("Plotting guide results")
  COMPOSE::plot_sigguides(results, cutoff = cutoff, sig.cutoff = sig.cutoff)
  for(i in 1:ncol(results$L2FC)){
    COMPOSE::plot_volcano(results, contrast = i, top = top.labels)
  }
  if(verbose)
    print("Calculating gene scores and plotting results")
  genescore <- COMPOSE::guides2genes(results$L2FC, neg.controls = neg.controls, essential.genes = essential.genes,
                            save.tables = save.table, split = split)

  COMPOSE::L2FC.violin(genescore, input = 'Gene Scores', variable = colnames(metadata)[meta.idnum[3]])

  intersections <- COMPOSE::plot_upset(genescore, list.return = list.return, cutoff = cutoff)

  if(save.plot)
    dev.off()
  if(pathway.analysis){
    if(verbose)
      print("Running Reactome pathway analyses")
    if(save.plot)
      pdf(paste0('pathway.analyses',Sys.Date(),'.pdf'), width = 12, height = 8)
    COMPOSE::reactomepa.enrich(genescore[!rownames(genescore) %in% essential.genes,],
                               genome.db = gene.db, showCategory = top.labels,
                               cutoff = cutoff, pvalue = sig.cutoff, save.table = save.table)
    for(i in 1:ncol(genescore)){
      COMPOSE::reactomepa.enrich(genescore[!rownames(genescore) %in% essential.genes,i,drop=FALSE],
                                 genome.db = gene.db, showCategory = top.labels,
                                 cutoff = cutoff, pvalue = sig.cutoff, save.table = save.table)
    }
    if(save.plot)
      dev.off()
  }
  results$genescore <- genescore
  return(results)
}


