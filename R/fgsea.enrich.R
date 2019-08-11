#' FGSEA utility
#'
#' fgsea.enrich utilizes the package 'fgsea' to test for enrichment of user-given genesets in your gene list
#'
#' @param CSS output CRISPR screen score matrix from create.CSS()
#' @param geneset input geneset
#' @param genome.db genome database from which to convert gene symbols to entrez ids
#' @param gene.ids type of gene ids in the geneset. Options are currently either 'symbol' or 'entrez'
#' @param p.value specified cutoff for pathway significance
#' @param minSize low cutoff for fgsea gene set size
#' @param maxSize high cutoff for fgsea gene set size
#' @param nperm number of permutations for fgsea analysis
#' @param top number of top pathways to display
#' @param collapsepath logical: collapse pathways using fgsea's collapsePathways function
#' @return object of fgsea output
#' @seealso \code{\link{fgsea}}
#' @export
#' @examples
#' fgsea.enrich <- reactomepa.enrich(CSS.data[,2, drop=FALSE])
#' fgsea.enrich <- reactomepa.enrich(CSS.data[,c(2:4)])
#' ...
#' @importFrom stats na.omit
#' @import fgsea

fgsea.enrich <- function(CSS, geneset, genome.db = Org.Hs.eg.db, gene.ids = "SYMBOL", p.value = 0.05, minSize = 15, maxSize = 500,
    nperm = 10000, top = 10, collapsepath = FALSE) {
  if (gene.ids != "SYMBOL") {
    CSS <- id.convert(CSS, genome.db = genome.db, gene.ids = gene.ids)
    # create vector for fgsea input (average CSS values)
    test <- rowMeans(CSS)
    names(test) <- rownames(CSS)
  }
  if (gene.ids == "SYMBOL") {
    # create vector for gage input (average CSS values)
    test <- rowMeans(CSS)
    names(test) <- rownames(CSS)
  }
    test <- test[order(test, decreasing = T)]
    fgseaRes <- fgsea::fgsea(pathways = geneset, stats = test, minSize = minSize, maxSize = maxSize, nperm = nperm)
    fgseaRes <- data.frame(fgseaRes[order(fgseaRes$pval), ])
    cat(sum(fgseaRes$padj < p.value), " significantly enriched gene sets", "\n")
    topPathwaysUp <- fgseaRes[fgseaRes$ES > 0, ]$pathway[1:top]
    topPathwaysDown <- fgseaRes[fgseaRes$ES < 0, ]$pathway[1:top]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    grid::grid.newpage()
    fgsea::plotGseaTable(geneset[topPathways], test, fgseaRes, gseaParam = 0.5)
    if (collapsepath) {
        collapsedPathways <- fgsea::collapsePathways(fgseaRes[order(pval)][padj < p.value], geneset, test)
        mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
        fgsea::plotGseaTable(geneset[mainPathways], test, fgseaRes, gseaParam = 0.5)
    }
    return(fgseaRes)
}


