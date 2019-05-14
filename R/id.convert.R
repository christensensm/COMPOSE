#' ID converter
#'
#' This function will take CSS with gene ids as rownames and convert them to a user-defined different ID using a user-given Annotation database
#'
#' @param CSS output CRISPR screen score matrix from create.CSS()
#' @param genome.db genome database from which to convert gene symbols to other ids
#' @param gene.ids type of gene ids in the geneset
#' @return object with converted rownames to the desired gene IDs
#' @seealso \code{\link[AnnotationDbi]{mapIds}}
#' @export
#' @examples
#' data(CSS.data)
#' library(org.Hs.eg.db)
#' CSS.newids <- id.convert(CSS.data, genome.db = genome.db, gene.ids = gene.ids)
#' ...
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit

id.convert <- function(CSS, genome.db = Org.Hs.eg.db, gene.ids = "SYMBOL"){
  CSS.new <- data.frame(CSS)
  CSS.new$NEWID = AnnotationDbi::mapIds(genome.db, keys = rownames(CSS.new),
                                        column = gene.ids, keytype = "SYMBOL", multiVals = "first")
  #Remove any symbols that return NA after conversion
  CSS.new <- stats::na.omit(CSS.new)
  rownames(CSS.new) <- CSS.new$NEWID
  CSS.new <- CSS.new[,-ncol(CSS.new), drop = FALSE]
  return(CSS.new)
}

