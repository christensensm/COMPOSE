#' gage Geneset Enrichment Test
#'
#' gage.enrich utilizes the package 'gage' to test for enrichment of user-input gene sets
#'
#' @param CSS output CRISPR screen score matrix from create.CSS()
#' @param geneset input geneset
#' @param genome.db genome database from which to convert gene symbols to entrez ids
#' @param gene.ids type of gene ids in the geneset. Options are currently either 'symbol' or 'entrez'
#' @param species species for pathview (if KEGG genesets are used)
#' @param database for visualizing pathways (use only if your database is KEGG)
#' @return list of output from gage
#' @seealso \code{\link{gage}, \link{pathview}}
#' @export
#' @examples
#' enriched.kegg <- gage.enrich(results$CSS, kegg.genesets, gene.ids = 'entrez',
#' species = 'hsa', kegg = T)
#' enriched.msigdb <- gage.enrich(results$CSS, c2.canonical, gene.ids = 'symbol')
#' ...
#' @importFrom stats na.omit
#' @importFrom gage gage


gage.enrich <- function(CSS, geneset, genome.db = org.Hs.eg.db, gene.ids = "SYMBOL", species = "hsa", database = NULL) {
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
    # same.dir signifies if changes are in a single direction (T) or both directions (F)
    FC.geneset <- gage::gage(test, gsets = geneset)
    FC.geneset$greater <- stats::na.omit(FC.geneset$greater)
    FC.geneset$less <- stats::na.omit(FC.geneset$less)
    sel.greater <- FC.geneset$greater[, "q.val"] < 0.1 & !is.na(FC.geneset$greater[, "q.val"])
    sel.less <- FC.geneset$less[, "q.val"] < 0.1 & !is.na(FC.geneset$less[, "q.val"])
    if (sum(sel.greater) > 0 | sum(sel.less) > 0) {
        path.ids <- rownames(FC.geneset$greater)[sel.greater]
        path.ids.l <- rownames(FC.geneset$less)[sel.less]
        if (!is.null(database)) {
            if (toupper(database) == "KEGG") {
                path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
                pv.out.list <- sapply(path.ids2[1:min(5, length(path.ids2))], function(pid) pathview::pathview(gene.data = test,
                  pathway.id = pid, species = species))
            }
        }
    }
    if (sum(sel.greater) == 0 & sum(sel.less) == 0)
        cat("Zero enriched genesets")
    return(FC.geneset)
}

