#' Create annotation object from .obo and .gaf files
#'
#' @param onto_path - full path to .obo file
#' @param annot_path - full path to .gaf file
#'
#' @return - list object contains annotation data ready for functional annnotation procedure,
#'           ontology relation data, GO terms full names, GO terms namespaces and metadata for annotation
#' @export
#'
#' @examples
#'     onto_path <- system.file("extdata", "goslim_plant.obo", package = "fsgor")
#'     annot_path <- system.file("extdata", "gene_association.tair", package = "fsgor")
#'     annot <- prepare_annotation(onto_path, annot_path)
#'
#'     # extract annotation calculated considering GO term ancestry
#'     annot$annot_data
#'
#'     # GO term ontology relations
#'     annot$relation_data
#'
#'     # GO terms real names
#'     annot$real_names
#'
#'     # GO terms namespaces
#'     annot$namespaces
#'
#'     # annotation metadata
#'     annot$annot_metadata
#'
prepare_annotaion <- function(onto_path, annot_path){
  outlist <- list()
  onto_data <- fsgor::read_obo(onto_path)
  annot_data <- fsgor::read_gaf(annot_path)
  annot_metadata <- annot_data$info

  # convert annotation from gaf format to list object
  # and create full annotation with ancestry according True Path Rule
  annot_data <- fsgor::convert_annotation_wide(annot_data$data)

  annot_data <- fsgor::create_annotation_with_ancestry(onto_data$relation_data, annot_data)

  outlist$annot_data <- annot_data
  outlist$relation_data <- onto_data$relation_data
  outlist$real_names <- onto_data$real_names
  outlist$namespaces <- onto_data$namespaces
  outlist$annotation_metadata <- annot_metadata
  return(outlist)
}
