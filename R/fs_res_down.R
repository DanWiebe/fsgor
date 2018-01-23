#' List contains dataframes with fold-specificty recognition results for example set of down-regualted genes
#'
#' @format A list with 2 entries
#' \describe{
#'     \item{fold-change-specific GO terms}{dataframe consitst of following columns: ids, namespace, name, padj, interval, genes}
#'     \item{not fold-change-specific GO terms}{dataframe consitst of following columns: ids, namespace, name, padj, interval, genes}
#' }
#'
#' @docType data
"fs_res_down"
