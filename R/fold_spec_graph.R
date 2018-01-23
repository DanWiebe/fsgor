#' Create igraph object for fold-specific terms
#'
#' @param fs_terms - vector contains fold-specific terms GO ids
#' @param relation_data - ontology relation data (prepare_annotation(...)$relation_data
#' @param reg - regulation type for DEG taken in to analysis ("up" or "down")
#'
#' @return - igraph object
#'
#' @examples
#' \dontrun{
#'    data(fs_res_up, annot_data)
#'    igraph_data_prep(fs_res_up$fs$ids, annot_data, "up")
#' }
igraph_data_prep <- function(fs_terms, relation_data, reg){

  terms_to_stay <- c(character(0))
  for (fs_term in fs_terms){
    terms_to_stay <- c(terms_to_stay,
                       fsgor::get_term_ancestors(fs_term,
                                                 relation_data$relation_data))
  }
  terms_to_stay <- unique(c(terms_to_stay, fs_terms))

  real_names <- relation_data$real_names
  relation_data <- relation_data$relation_data[terms_to_stay]

  from <- c(character(0))
  to <- c(character(0))

  for (name in names(relation_data)){
    for (value in relation_data[[name]]){
      to <- c(to, name)
      from <- c(from, value)
    }
  }

  actors <- data.frame(name = unique(union(to, from)))
  relations <- data.frame(from = from, to = to)
  graph <- igraph::graph_from_data_frame(relations,
                                         directed = TRUE, vertices = actors)
  igraph::V(graph)$shape <- "circle"
  igraph::V(graph)$label.color <- "black"
  igraph::V(graph)$label.degree <- pi / 2
  igraph::V(graph)$label.dist <- 0.45
  igraph::V(graph)$size <- 5

  if (reg == "up") {
    igraph::V(graph)[fs_terms]$color <- "yellow2"
  } else if (reg == "down"){
    igraph::V(graph)[fs_terms]$color <- "skyblue2"
  }
  igraph::V(graph)[setdiff(terms_to_stay, fs_terms)]$color <- "gray"
  for (term in terms_to_stay){
    real_name_vec <- unlist(strsplit(real_names[[term]], "\\s"))
    wrapped__name <- paste(real_name_vec, collapse = "\n")
    igraph::V(graph)[[term]]$name <- paste(term, wrapped__name, sep = "\n")
  }
  return(graph)
}

#' Extract all ancestors for distinct GO term
#'
#' @param term_vec - vector of terms
#' @param relation_data - ontology relation data (prepare_annotation(...)$relation_data
#'
#' @return
#'
#' @examples
#' \dontrun{
#'     onto_data <- fsgor::read_obo("~/go.obo")
#'     get_term_ancestors(term_vec, onto_data$relation_data)
#' }
get_term_ancestors <- function(term_vec, relation_data){
  output_vec <- c(character(0))
  while (length(term_vec) != 0){
    temp_vec <- c(character(0))
    for (term in term_vec){
      temp_vec <- c(temp_vec, relation_data[[term]])
    }
    term_vec <- unique(temp_vec)
    output_vec <- c(output_vec, term_vec)
  }
  return(unique(output_vec))
}

#' Create igraph object for fold-specificity recognition function output (recognize_fs_terms())
#'
#' @param fs_res - list contains fold-specificity recognition results
#' @param annot_data - prepare_annotation() function output
#' @param reg - regulation type for DEG taken in to analysis ("up" or "down")
#' @param scope - GO term namespace ("biological_process", "cellular_component", "molecular_function"),
#'                NULL by default (GO terms from all namespaces will be plotted)
#'
#' @return - igraph object with sugiyama layout
#' @export
#'
#' @examples
#'      data(fs_res_up, annot_data, package = "fsgor")
#'      fold_spec_graph(fs_res_up, annot_data, "up", scope = "molecular_function")
fold_spec_graph <- function(fs_res, annot_data, reg, scope = NULL){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  fs_terms <- fs_res$fs$ids
  if (!is.null(scope)) {
    fs_terms <- (fs_res$fs[grep(scope, fs_res$fs$namespace, fixed = TRUE), ])$ids
  } else {
    fs_terms <- fs_res$fs$ids
  }


  g <- fsgor::igraph_data_prep(fs_terms = fs_terms,
                               relation_data = annot_data, reg = reg)
  g$layout <- igraph::layout_with_sugiyama(g)$layout
  return(g)
}
