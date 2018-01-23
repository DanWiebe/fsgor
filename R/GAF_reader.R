#' GAF file reader
#'
#' @param gafpath - full path to GAF file
#'
#' @return - list contains information about file and data accessible by $info and $data keys correspondingly
#'
#' @examples
#'     gaf_path <- system.file("extdata", "gene_association.tair", package = "fsgor")
#'     read_gaf(gaf_path)
read_gaf <- function(gafpath){
  lines <- readLines(gafpath)

  # separate header and data parts of file
  inds <- grep("!.*", lines)
  header <- lines[inds]
  lines <- lines[-inds]

  # split lines and put them in data frame with specified col names
  df <- matrix(character(0), ncol = 17, nrow = length(lines))
  df <- t(sapply(lines, function(x) strsplit(x, "\t")[[1]]))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  colnames(df) <- c("DB", "DB_Object_ID", "DB_Object_Symbol",
                    "Qualifier", "GO_ID", "DB:Reference",
                    "Evidence_Code", "With_(or)_From", "Aspect",
                    "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
                    "Taxon", "Date", "Assigned_By",
                    "Annotation_Extension", "Gene_Product_Form_ID")
  outlist <- list(info = header, data = df)
  return(outlist)
}
