#' Convert GAF format type annotation to list contains GO term id's as keys and Gene ID's as values
#'
#' @param data - GAF annotation table (output of read.gaf function)
#' @param id_pattern - regex for Gene identifier (see fsgor::patterns_enum$Arabidopsis_thaliana for Arabidopsis AGI code example pattern)
#'
#' @return - list with GO term id's as keys and Gene ID's as values
#'
#' @examples
#' \dontrun{
#'     convert_annotation_wide(annot_data, fsgor::patterns_enum$Arabidopsis_thaliana)
#' }
convert_annotation_wide <-
  function(data, id_pattern = NULL) {
    # extract GO id's and Gene id's annotated to them
    annotdf <-
      data.frame(
        "GOID" = data$GO_ID,
        "GeneID" = data$DB_Object_Name,
        stringsAsFactors = FALSE
      )
    # leave Gene ID's only of appropriate format if id_pattern is not NULL
    if (!is.null(id_pattern)) {
      annotdf <- annotdf[grep(id_pattern, annotdf[, 2]), ]
    }
    # convert data frame to list with GO term id's as keys
    # and corresponding gene id's as values
    annotdf[["GOID"]] <- as.factor(annotdf[["GOID"]])
    outlist <- split(annotdf, annotdf[["GOID"]])
    outlist <- lapply(outlist, function(x)
      x[, -1])
    return(outlist)
  }
