#' Extraction of gene names
#'
#' Extraction of gene names from rownames that contain single guide RNA IDs
#'
#' @param data input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param neg.controls input dataframe containing sample names and other identifiers or data
#' @param split string identifying cell line on which to filter counts
#' @return new data with added column of gene names
#' @export
#' @examples
#' new.L2FC <- extract.genefromguide(L2FC, neg.controls, split = '_')


extract.genefromguide <- function(data, neg.controls, split = "_") {
    # Pull out gene names from guide ids and add as a column to data
    data <- data[order(rownames(data)), ]
    gRNAnames <- rownames(data)
    gRNA.neg <- gRNAnames[gRNAnames %in% neg.controls]
    gRNA.noneg <- gRNAnames[!gRNAnames %in% neg.controls]
    gRNAgenes <- sapply(as.character(gRNA.noneg), function(x) strsplit(x, split = split)[[1]][1])
    gRNAgenes <- c(gRNAgenes, gRNA.neg)
    names(gRNAgenes) = NULL
    gRNAgenes = gRNAgenes[order(gRNAgenes)]
    new.data <- data.frame(data)
    new.data$gene.name <- gRNAgenes
    return(new.data)
}
