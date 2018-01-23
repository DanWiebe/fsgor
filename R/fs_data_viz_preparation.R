#' Create rectangle plot coordinates for each Go term
#'
#' @param intervals - vector of intervals in character representation (e.g. "1-6")
#'
#' @return - matrix contains interval start coordinates in the first column and interval end coordinates
#'           in the second row (note: the start coordinate is deducted by 1 in order to make the appropriate layout for rectangles)
#'
#' @examples
#' \dontrun{
#'     create_intervals_matrix(intervals)
#' }
create_intervals_matrix <- function(intervals) {
  output_matrix <- matrix(integer(0), ncol = 2)
  # convert character representation of coordinates to numeric
  for (i in 1:length(intervals)) {
    x <- intervals[i]
    if (grepl("^([1-9][0-9]*)-([1-9][0-9]*)$", x) == TRUE) {
      splitarr <- unlist(strsplit(x, "-", fixed = TRUE))
      start <- as.numeric(splitarr[1]) - 1
      end <- as.numeric(splitarr[2])
      if (start > end)
        stop (paste("start position of interval is greater than end position:", x))
      output_matrix <- rbind(output_matrix, c(start, end))
    } else if (grepl("^([1-9][0-9]*)$", x) == TRUE) {
      start <- as.numeric(x) - 1
      end <- as.numeric(x)
      output_matrix <- rbind(output_matrix, c(start, end))
    } else {
      start <- as.numeric(x)
      end <- as.numeric(x)
      output_matrix <- rbind(output_matrix, c(start, end))
    }
  }
  return(output_matrix)
}

#' Data parser for fold-specificity rectangle plot
#'
#' Create input data for rectangle plot (fold_spec_chart function) using recognize_fs_terms function output as input
#'
#' @param fs_res_up dataframe contains fold-specificty recognition data (recognize_fs_terms function output) for up regulation
#' @param fs_res_down dataframe contains fold-specificty recognition data (recognize_fs_terms function output) for down regulation
#' @param scope namespace of terms for plot input data (NULL by default)
#'
#' @return input dataframe for create.fold.spec.chart function
#'
#' @examples
#'     data(fs_res_up, fs_res_down, package = "fsgor")
#'     fold_spec_chart_data(fs_res_up, fs_res_down, scope = "molecular_function")
fold_spec_chart_data <-
  function(fs_res_up, fs_res_down, scope = NULL) {

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    # leave data only for definite namespace if scope is not NULL
    if (!is.null(scope)) {
      fs_res_up <-
        fs_res_up[grep(scope, fs_res_up$namespace, fixed = TRUE), ]
      fs_res_down <-
        fs_res_down[grep(scope, fs_res_down$namespace, fixed = TRUE), ]
    }

    # add column with regulation type
    fs_res_up$reg <- rep("up", length(fs_res_up[, 1]))
    fs_res_down$reg <- rep("down", length(fs_res_down[, 1]))

    # combine up and down dataframes and
    # spread the resulting one by regulation type
    fs_res_combined <- rbind(fs_res_up, fs_res_down)
    fs_res_combined <- fs_res_combined[, c(-1, -2, -4)]
    fs_res_combined_spreaded <-
      tidyr::spread(fs_res_combined, "reg", "interval", fill = "0")

    # add zero's if there is no data for certain regulation type
    if (is.null(fs_res_combined_spreaded$up)) {
      fs_res_combined_spreaded$up <-
        rep("0", length(fs_res_combined_spreaded$name))
    }

    if (is.null(fs_res_combined_spreaded$down)) {
      fs_res_combined_spreaded$down <-
        rep("0", length(fs_res_combined_spreaded$name))
    }

    # create matrices with intervals coordinates
    in_mat_up <- fsgor::create_intervals_matrix(fs_res_combined_spreaded$up)
    in_mat_down <- fsgor::create_intervals_matrix(fs_res_combined_spreaded$down)

    # add start - end coordinates to resulting dataframe
    fs_res_combined_spreaded$start_up <- in_mat_up[, 1]
    fs_res_combined_spreaded$end_up <- in_mat_up[, 2]
    fs_res_combined_spreaded$start_down <- in_mat_down[, 1]
    fs_res_combined_spreaded$end_down <- in_mat_down[, 2]
    fs_res_combined_spreaded <- fs_res_combined_spreaded[, c(-2, -3, -4)]
    fs_res_combined_spreaded <-
      fs_res_combined_spreaded[order(fs_res_combined_spreaded[, 3]), ]

    return(fs_res_combined_spreaded)
}
