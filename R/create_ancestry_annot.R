#' Annotation with ancestry creator
#'
#' Annotaion for terms taking into account True Path Rule
#'
#' @param rellist list contains child terms as keys and nearest parents as values
#' @param annotlist list contains terms as keys and direct annotations as values
#'
#' @return list contains terms as keys and genes annotated according to True Path Rule as values
#'
#' @examples
#' \dontrun{
#'     onto_data <- fsgor::read_obo("~/go.obo")
#'     annot_data <- fsgor::read_gaf("~/goa_human.gaf")
#'     annot_data <- fsgor::convert_annotation_wide(annot_data$data)
#'     create_annotation_with_ancestry(onto_data$relation_data, annot_data)
#' }
#'
create_annotation_with_ancestry <- function(rellist, annotlist) {
  # extract terms that are leaves and roots in GO DAG to corresponding vectors
  leaves <- fsgor::find_boundaries(rellist)$leaves
  roots <- fsgor::find_boundaries(rellist)$roots

  # delete terms that are not presented in relation file
  # (it must be a list entry)
  terms_not_in_rel <-  setdiff(names(annotlist), names(rellist))
  annotlist <-
    annotlist[-match(setdiff(terms_not_in_rel, roots), names(annotlist))]

  # traversing graph in down-up manner starting from leaves to the roots
  while (length(leaves) != 0) {
    #temporary vector for storing parental terms for next iteration of loop
    tempvec <- c()
    for (i in 1:length(leaves)) {
      if (!leaves[i] %in% roots) {
        # put all nearest parents for current term in parents vector
        # then iterate through parents vector
        # and assign them child term annotation
        parents <- rellist[[leaves[i]]]
        for (j in 1:length(parents)) {
          annotlist[[parents[j]]] <-
            unique(c(annotlist[[parents[j]]], annotlist[[leaves[i]]]))
        }

        # append parents of current term to temporary vector
        tempvec <- unique(c(tempvec, parents))
      }
    }

    # reassign leaves vector with parents for next iteration of loop
    leaves <- tempvec
  }
  return(annotlist)
}

#' DAG roots and leaves extractor
#'
#' @param rellist list contains child terms as keys and nearest parents as values
#'
#' @return list contains root and leaf terms
#'
#' @examples
#' \dontrun{
#'     onto_data <- fsgor::read_obo("~/go.obo")
#'     find_boundaries(onto_data$relation_data)
#' }
find_boundaries <- function(rellist) {
  outlist <- list()

  # extract child terms and parental terms
  children <- names(rellist)
  parents <- unlist(rellist)

  # define terms that are not  presented in the parental set as leaves
  outlist$leaves <- setdiff(children, parents)

  # define terms that are not presented in the children set as roots
  outlist$roots <- setdiff(parents, children)

  return(outlist)
}
