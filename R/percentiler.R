#' Gene set quantile separator
#'
#' Converts the set of genes in to quantiles according log value of fold change and
#' produces all possible combinations of resulting gene subsets.
#' Input dataset must be dataframe with gene ID's in the first row and corresponding
#' fold change values in the second row. Note: dataframe must contain only genes with
#' common type of regulation (up or down)
#'
#' @param data dataframe contains initial set of genes gene ID's in the first row
#'             and corresponding fold change values in the second row
#' @param n number of quantiles
#'
#' @return list of gene sets for each quatile and all combinations
#' @export
#'
#' @examples
#'     data(up_genes, package = "fsgor")
#'     gene_sets <- div_genes_to_quantiles(up_genes, 6)
#'
div_genes_to_quantiles <- function(data, n) {
  data <- data[order(data[, 2]), ]
  folds <- data[, 2]
  genes <- data[, 1]
  percentiles <- c(1:n) / n
  borders <- stats::quantile(folds, percentiles)
  negcheck <- sum(folds < 0)
  singleintervals <- list()
  if (negcheck == 0) {
    borders <- c(folds[1], borders)
    for (i in 1:(n - 1)) {
      singleintervals[[toString(i)]] <-
        genes[folds >= borders[i] & folds < borders[i + 1]]
    }
    singleintervals[[toString(n)]] <-
      genes[folds >= borders[length(borders) - 1]]
  } else if (negcheck == length(folds)) {
    folds <- rev(folds)
    genes <- rev(genes)
    borders <- rev(borders)
    for (i in 1:(n - 1)) {
      singleintervals[[toString(i)]] <-
        genes[folds <= borders[i] & folds > borders[i + 1]]
    }
    singleintervals[[toString(n)]] <-
      genes[folds <= borders[length(borders)]]
  } else {
    stop(
      "input data contains both up and down regulated genes, please analyze them separately",
      call. = FALSE
    )
  }
  crossintervals <- list()
  for (i in 1:(n - 1)) {
    finalvec <- singleintervals[[toString(i)]]
    for (j in (i + 1):n) {
      name <- paste(i, j, sep = "-")
      finalvec <- c(finalvec, singleintervals[[toString(j)]])
      crossintervals[[name]] <- unlist(finalvec)
    }
  }
  return(append(singleintervals, crossintervals))
}
