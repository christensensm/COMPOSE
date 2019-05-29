#' Fold-change calculation for CRISPR screen
#'
#' Calculates log2 fold-changes between samples of a specified variable using input of a count matrix and metadata dataframe
#'
#' @param countsTable input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param metadata input dataframe containing sample names and other identifiers or data
#' @param meta.idnum vector containing the column numbers in metadata that represent (1) Cell line, (2) replicates, and (3) condition
#' @param include.batch logical to include replicates in the model
#' @param plot.volcano logical - do you want to show the volcano plot for each cell type
#' @param p.cutoff numeric - specified p-value cut-off
#' @param plot.MA logical - do you want to show the MA plot for each cell type
#' @param save logical - do you want to save the fold-change table to csv
#' @return matrix containing Log2 fold-changes for each comparison
#' @export
#' @examples
#' L2FC <- calc.DESeq2.L2FC(countsTable, design.table, plot.MA = T, save = T)
#' ...
#' @importFrom utils combn write.csv
#' @import BiocParallel
#' @import parallel
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink plotMA resultsNames
#' @importFrom stats density

# Create function to calculate Log2 fold change for each 'sample' between all 'condition' factors
calc.DESeq2.L2FC <- function(countsTable, metadata, meta.idnum = NULL, include.batch = F, plot.volcano = F, p.cutoff = 0.05,
    plot.MA = F, save = F) {
    if (ncol(countsTable) != nrow(metadata))
        stop("Differing number of samples in count matrix and metadata table")
    if (!all(colnames(countsTable) == metadata[, 1]))
        stop("Please make sure your count matrix column names are the same as your metadata sample IDs (first column)")
  design.table <- data.frame(metadata[, meta.idnum])
  colnames(design.table) <- c("sample", "replicates", "condition")
  design.table[, 3] <- factor(design.table[, 3], levels = unique(design.table[, 3]))  #condition factor
  ncontrasts <- ncol(utils::combn(unique(design.table[, 3]), 2)) * length(unique(design.table[, 1]))
  L2FC <- matrix(0, 0, ncol = ncontrasts, nrow = nrow(countsTable))
  colnames(L2FC) <- character(length = ncol(L2FC))
  padj <- matrix(0, 0, ncol = ncontrasts, nrow = nrow(countsTable))
  colnames(padj) <- character(length = ncol(padj))
  lfcSE <- matrix(0, 0, ncol = ncontrasts, nrow = nrow(countsTable))
  colnames(lfcSE) <- character(length = ncol(lfcSE))
  current.col = 1
  graphics::par(mfrow = c(2, 2))
  for (i in 1:length(unique(design.table[, 1]))) {
    # go through each cell sample
    print(paste0("Calculating DESeq2 fold-changes for ", unique(design.table[, 1])[i]))
    # create counts/design.table table of cell sample i
    temp.counts <- countsTable[, which(design.table[, 1] == unique(design.table[, 1])[i])]
    temp.design <- design.table[which(design.table[, 1] == unique(design.table[, 1])[i]), ]
    condition.fac <- factor(temp.design[, 3], levels = unique(temp.design[, 3]))  #condition factor
    contrast.cond <- combn(unique(condition.fac), 2)
    if (include.batch)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = temp.counts,
                                            colData = temp.design, design = ~replicates + condition)
    if (!include.batch)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = temp.counts, colData = temp.design, design = ~condition)
    workers = parallel::detectCores() - 1
    paraparam = BiocParallel::MulticoreParam(workers = workers, exportglobals = FALSE)
    dds <- DESeq2::DESeq(dds, parallel = T, BPPARAM = paraparam, quiet = T)
    for (j in 1:ncol(contrast.cond)) {
      res <- DESeq2::results(dds, contrast = c("condition", as.character(contrast.cond[2, j]),
                                               as.character(contrast.cond[1, j])))
      resOrdered <- as.data.frame(res[order(rownames(res)), ])
      L2FC[, current.col] <- resOrdered$log2FoldChange
      L2FC[, current.col] <- scale(L2FC[,current.col])
      padj[, current.col] <- resOrdered$padj
      padj[, current.col][is.na(resOrdered$padj)] = 1
      lfcSE[, current.col] <- resOrdered$lfcSE
      naming <- paste0(unique(design.table[, 1])[i], "_",
                       as.character(contrast.cond[2, j]), "_vs_",
                       as.character(contrast.cond[1,j]))
      if (save){
        reswrite <- resOrdered[order(resOrdered$log2FoldChange), ]
        utils::write.csv(as.data.frame(reswrite),
                         file = paste0(naming, "_results.csv"))
      }
      colnames(L2FC)[current.col] <- colnames(padj)[current.col] <- colnames(lfcSE)[current.col] <- naming
      current.col <- current.col + 1
    }
  }
  graphics::par(mfrow = c(1, 1))
  rownames(L2FC) <- rownames(padj) <- rownames(lfcSE) <- rownames(resOrdered)

  results <- list(countsTable, L2FC, padj, lfcSE)
  names(results) <- c('counts','L2FC','padj','lfcSE')
  return(results)
}
