#' Functional annotation test
#'
#' Conducts functional annotation using Fisher exact test
#'
#' @param annotlist - list contains terms as keys and direct annotations as values
#' @param realnames - GO terms full names
#' @param namespaces - GO terms namespaces
#' @param inputgeneset - vector of query genes
#' @param gene_amount_border - minimal number of genes annotated to a term (1 by default)
#' @param fisher_alternative - fisher exact test alternative,
#'                           possible values: "greater", "less", "two.sided" ("greater by default)
#' @param p_adjust_method - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'
#' @return - dataframe contains functional annotaion for gene set
#'
#' @examples
#'     # data - dataframe contains gene IDs in first column and fold-change values in second column
#'     data(up_genes, annot_data, package = "fsgor")
#'     gene_sets <- div_genes_to_quantiles(up_genes, 6)
#'     func_annot_test_internal(annot_data$annot_data, annot_data$real_names, annot_data$namespaces, gene_sets$'1')
func_annot_test_internal <- function(annotlist,
                                     realnames,
                                     namespaces,
                                     inputgeneset,
                                     gene_amount_border = 1,
                                     fisher_alternative = "greater",
                                     p_adjust_method = "BH") {

  # extract background genes set from annotation list
  # and calculate the amount
  bggenes <- unique(unlist(annotlist))
  bgtot <- length(bggenes)

  # delete genes from query genes set that are not presented
  # in the background and calculate the amount
  inputgeneset <- intersect(inputgeneset, bggenes)
  qtot <- length(inputgeneset)

  annotlist_names <- names(annotlist)

  qit_list <-
    lapply(annotlist, function(x)
      intersect(inputgeneset, x))

  names(qit_list) <- annotlist_names

  qit_ints <- lengths(qit_list)

  qit_list <- qit_list[qit_ints >= gene_amount_border]

  qit_ints <- lengths(qit_list)

  annotlist <- annotlist[match(names(qit_list), names(annotlist))]

  annotlist_names <- names(annotlist)

  bg_ints <- lengths(annotlist)

  values_matrix <-
    cbind(qit_ints, rep(qtot, length(qit_ints)),
          bg_ints, rep(bgtot, length(qit_ints)))

  fisher_matrix <-
    cbind(
      values_matrix[, 1],
      values_matrix[, 2] - values_matrix[, 1],
      values_matrix[, 3] - values_matrix[, 1],
      values_matrix[, 4] - values_matrix[, 2]
      - values_matrix[, 3] + values_matrix[, 1]
    )

  pvalues <-
    apply(fisher_matrix, 1, function(x)
      stats::fisher.test(matrix(x, nrow = 2), fisher_alternative)$p.value)

  real_names <-
    unlist(realnames[match(names(annotlist), names(realnames))])

  name_spaces <-
    unlist(namespaces[match(names(annotlist), names(namespaces))])

  genes <-
    unlist(lapply(qit_list, function(x)
      paste(x, collapse = "/")))

  padj <- stats::p.adjust(pvalues, method = p_adjust_method)
  outdf <- data.frame(
    "GO_id" = annotlist_names,
    "namespace" = name_spaces,
    "name" = real_names,
    "qit" = values_matrix[, 1],
    "qtot" = values_matrix[, 2],
    "bgit" = values_matrix[, 3],
    "bgtot" = values_matrix[, 4],
    "pval" = pvalues,
    "padj" = padj,
    "qit_genes" = genes,
    stringsAsFactors = FALSE
  )

  # sort output df by p-values in ascending order
  outdf <- outdf[order(outdf$pval), ]
  # add row numbering to output data frame
  rownames(outdf) <- seq(nrow(outdf))
  return(outdf)
}
