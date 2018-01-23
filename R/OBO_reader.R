#' OBO file reader
#'
#' Parses .obo format file (open biomedical ontologies), returns list contains GO term ontology relation data,
#' GO terms real names, GO terms namespaces
#'
#'
#' @param obopath path to the .obo format file
#'
#' @return lists that contain GO terms ancestry data, names for GO terms, and corresponding namespace for each GO term
#'
#' @examples
#'     go_path <- system.file("extdata", "goslim_plant.obo", package = "fsgor")
#'     read_obo(go_path)
read_obo <- function(obopath) {
  # read all lines from file
  lines <- readLines(obopath)
  # empty lists for data and id vector for storing GO ids
  outdata <- list()
  datalist <- list()
  nameslist <- list()
  namespacelist <- list()
  id_vector <- c(character(0))
  # counters: i - for iterating trough lines,
  #           ind - for creating indices for list entries
  i <- 1
  list_ind <- 1
  while (i != length(lines)) {
    id <- ""
    name <- ""
    namespace <- ""
    # flag that assigned with "true" value if current GO term marked as obsolete
    is_obsolete <- "false"
    # vector for storing parental GO terms
    parents <- c(character(0))
    # process piece of lines that contain informaton about distinct GO term
    if (lines[i] == "[Term]") {
      i <- i + 1
      # extract GO term id
      if (startsWith(lines[i], "id:")) {
        id <- lines[i]
        id <- gsub("id:", "", id)
        id <- gsub("^\\s+|\\s+$", "", id)
        i <- i + 1
      }
      # extract GO term real name
      if (startsWith(lines[i], "name:")) {
        name <- gsub("name:", "", lines[i])
        name <- gsub("^\\s+|\\s+$", "", name)
        i <- i + 1
      }
      # extract GO term corresponding namespace
      if (startsWith(lines[i], "namespace:")) {
        namespace <- gsub("namespace:", "", lines[i])
        namespace <- gsub("^\\s+|\\s+$", "", namespace)
        i <- i + 1
      }
      # process the rest part of lines that refer to
      # current GO term and extract parental GO term ids
      # by "is_a" relation keyword or "true" value for is_obsolete flag
      while (lines[i] != "") {
        if (startsWith(lines[i], "is_a:")) {
          parent <- gsub("is_a:|!.*", "", lines[i])
          parent <- gsub("^\\s+|\\s+$", "", parent)
          parents <- c(parents, parent)
        } else if (startsWith(lines[i], "is_obsolete:")) {
          is_obsolete <- gsub("is_obsolete:", "", lines[i])
          is_obsolete <- gsub("^\\s+|\\s+$", "", is_obsolete)
        }
        i <- i + 1
      }
      # if GO term is not obsolete fill in the lists with the appropriate values
      if (is_obsolete == "false") {
        id_vector[[list_ind]] <- id
        datalist[[list_ind]] <- parents
        nameslist[[list_ind]] <- name
        namespacelist[[list_ind]] <- namespace
        list_ind <- list_ind + 1
      }
    }
    i <- i + 1
  }
  # assign GO ids as names to lists
  names(datalist) <- id_vector
  names(nameslist) <- id_vector
  names(namespacelist) <- id_vector
  # remove entries without values from relation data list
  datalist <- datalist[lapply(datalist, length) > 0]
  outdata$relation_data <- datalist
  outdata$real_names <- nameslist
  outdata$namespaces <- namespacelist
  return(outdata)
}
