#' Reactome pathway analysis
#'
#' reactome.enrich utilizes the package 'ReactomePA' to test for enrichment of reactome pathways in your gene list
#'
#' @param CSS output CRISPR screen score matrix from create.CSS()
#' @param genome.db genome database from which to convert gene symbols to entrez ids
#' @param pvalue specified cutoff for pathway significance (default = 0.05)
#' @param combine logical: do you want to combine the columns in the data? (default=FALSE)
#' @param combine.by vector separating the columns of your data into groups for comparing multiple groups
#' @param metric defines whether to look for enrichment based on a gene list of differentially expressed genes defined by a 'cutoff' (default) or by 'gsea' (one sample or fully combined samples only)
#' @param cutoff numeric cutoff to define differential expression (default = 1)
#' @param save logical: saves plots to pdf
#' @param showCategory number of categories to display
#' @param font.size font size of labels in dotplot and barplot
#' @return list of (1) output from ReactomePA enrichment tests and (2) readable output
#' @seealso \code{\link{ReactomePA}, \link{clusterProfiler}}
#' @export
#' @examples
#' reactome.enrich <- reactomepa.enrich(CSS.data[,1, drop=FALSE])
#' reactome.enrich <- reactomepa.enrich(CSS.data[,c(1:5)])
#' reactome.enrich <- reactomepa.enrich(CSS.data[,c(1:5)], combine = T)
#' reactome.enrich <- reactomepa.enrich(CSS.data[,1, drop = FALSE], metric = 'gsea')
#' ...
#' @importFrom ReactomePA enrichPathway gsePathway
#' @importFrom clusterProfiler compareCluster
#' @import grDevices
#' @importFrom stats na.omit
#' @importFrom enrichplot dotplot emapplot cnetplot

reactomepa.enrich <- function(CSS, genome.db = org.Hs.eg.db, pvalue = 0.05, combine = FALSE, combine.by = NULL, metric = "cutoff",
    cutoff = 1, save = FALSE, showCategory = 15, font.size = 8) {
  if (save)
    grDevices::pdf("ReactomePA_output.pdf", height = 8, width = 12)
  CSS <- id.convert(CSS, genome.db = genome.db, gene.ids = "ENTREZID")
  if (metric == "cutoff") {
        if (combine) {
            if (is.null(combine.by)) {
                # one list at a time based on gene ids
                test <- as.numeric(rowMeans(CSS[, 1:(ncol(CSS) - 1)]))
                names(test) <- rownames(CSS)
                test <- test[abs(test) >= cutoff]
                x <- ReactomePA::enrichPathway(names(test), pvalueCutoff = pvalue, readable = T)
                print(barplot(x, showCategory = showCategory, font.size = font.size))
                print(enrichplot::dotplot(x, showCategory = showCategory, font.size = font.size))
                print(enrichplot::emapplot(x))
                print(enrichplot::cnetplot(x, categorySize = "pvalue", foldChange = test) +
                        scale_color_gradient2(low = 'red', high = 'green', mid = 'white',
                                              limits = c(-max(abs(test)),max(abs(test)))))
            }
            if (!is.null(combine.by)) {
                test <- list()
                for (i in 1:length(unique(combine.by))) {
                  tempcombine <- which(combine.by == unique(combine.by)[i])
                  test[[i]] <- rownames(CSS)[abs(rowMeans(CSS[, tempcombine][, -ncol(CSS)])) >= cutoff]
                }
                names(test) <- unique(combine.by)
                x <- clusterProfiler::compareCluster(test, fun = "ReactomePA::enrichPathway")
                print(enrichplot::dotplot(x, showCategory = showCategory, font.size = font.size))
            }
        }
        if (!combine) {
            if (ncol(CSS) > 1) {
                test <- list()
                for (i in 1:ncol(CSS)) {
                  test[[i]] <- rownames(CSS)[abs(CSS[, i]) >= cutoff]
                }
                names(test) <- colnames(CSS)
                x <- clusterProfiler::compareCluster(test, fun = "ReactomePA::enrichPathway")
                print(enrichplot::dotplot(x, showCategory = showCategory, font.size = font.size) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
            }
            if (ncol(CSS) == 1) {
                test <- as.numeric(CSS[, 1])
                names(test) <- rownames(CSS)
                test <- test[abs(test) >= cutoff]
                x <- ReactomePA::enrichPathway(names(test), pvalueCutoff = pvalue, readable = T)
                print(barplot(x, showCategory = showCategory, font.size = font.size))
                print(enrichplot::dotplot(x, showCategory = showCategory, font.size = font.size))
                print(enrichplot::emapplot(x))
                print(enrichplot::cnetplot(x, categorySize = "pvalue", foldChange = test) +
                        scale_color_gradient2(low = 'red', high = 'green', mid = 'white',
                                              limits = c(-max(abs(test)),max(abs(test)))))
            }
        }
    }
  if (metric == "gsea") {
    test <- as.numeric(rowMeans(CSS))
    names(test) <- rownames(CSS)
    test <- test[test!=0]
    test <- test[order(test, decreasing = T)]
    x <- ReactomePA::gsePathway(test, nPerm = 10000, pvalueCutoff = pvalue, pAdjustMethod = "BH", verbose = FALSE)
    res <- as.data.frame(x)
    if (nrow(res) > 0)
      print(enrichplot::emapplot(x, color = "pvalue"))
  }
  res <- list(raw = CSS, output = x, readable = as.data.frame(x))
  if (save)
    grDevices::dev.off()
  return(res)
}

