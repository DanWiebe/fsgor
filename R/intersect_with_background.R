#' Intersect annotation with background set of genes
#'
#' Intersect list contains terms as keys and direct annotations as values
#' with vector contains background set of genes
#'
#' @param annot_data list contains terms as keys and direct annotations as values
#' @param background_genes vector contains background set of genes
#'
#' @return list contains terms as keys and direct annotations as values with genes only from background
#'
#' @examples
#'     data(annot_data, bg_genes, package = "fsgor")
#'     intersect_with_background(annot$annot_data, bg_genes)
intersect_with_background <-
  function(annot_data, background_genes) {
    term_ids <- names(annot_data)

    # intersect annotations for each term with background
    annot_data <-
      lapply(annot_data, function(x)
        intersect(x, background_genes))

    names(annot_data) <- term_ids

    # delete empty entries
    annot_data <- annot_data[lapply(annot_data, length) > 0]
    return(annot_data)
  }
