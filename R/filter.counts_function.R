#' Filtration of gRNAs with Low Counts
#'
#' Filtration and visualization of counts
#'
#' @param counts input matrix containing normalized gRNA counts with gRNA ids as row names
#' @param thresh numerical cutoff for low counts
#' @return new matrix of filtered counts
#' @export
#' @examples
#' y <- matrix(rnorm(100*9, mean = 10, sd = 1),100,9)
#' y[seq(1,100,10),] <- y[seq(1,100,10),]/20
#' counts.filtered <- filterCounts(y, thresh = 1, minsample.ids = 3)
#' ...
#' @importFrom edgeR cpm

# Filter counts and visualize filter via L2FC vs average counts requires metadata column 'sample' that identifies cell
# line Display normcounts vs average day 0 or control sample.ids to pick a cutoff point for excluding low readcount
# sgRNAs from all sample.ids.
filterCounts = function(counts, thresh = 1, minsample.ids = 2) {
    new.counts <- edgeR::cpm(counts)
    keep <- rowSums(new.counts > thresh) >= minsample.ids
    counts.filtered <- counts[keep, ]
    return(counts.filtered)
}
