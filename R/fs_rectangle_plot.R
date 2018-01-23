#' fold specificity rectangle plot
#'
#' Create rectangle plot for fold specificity data representation
#'
#' @param interval_labels vector of user defined names for non-overlaping intervals (e.g. c("weak response",...,"strong response")), see fsgor::six_perc_labs for example
#' @param x_text_size size of text for x axis labels
#' @param fs_res_up dataframe with fold-specific GO terms data (recognize_fs_terms() function output) for up regulated genes
#' @param fs_res_down dataframe with fold-specific GO terms data (recognize_fs_terms() function output) for down regulated genes
#' @param scope GO term namespace ("biological_process", "molecular_function", "cellular_component")
#'
#' @return fold specificity rectangle plot as ggplot object
#' @export
#' @importFrom ggplot2 ggplot geom_rect scale_x_continuous scale_y_continuous theme geom_hline geom_text coord_flip aes element_blank element_text
#'
#' @examples
#'     data(fs_res_up, fs_res_down, package = "fsgor")
#'     fold_spec_chart(fs_res_up, fs_res_down, interval_labels = fsgor::six_perc_labs)
fold_spec_chart <- function(fs_res_up,
                            fs_res_down,
                            scope = NULL,
                            interval_labels,
                            x_text_size = 10) {

  data <- fsgor::fold_spec_chart_data(fs_res_up$fs, fs_res_down$fs, scope)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # extract y-axis labels
  labels <- as.vector(data[, 1])
  # extract reactangle plot coordinates
  df <- data.frame(c(data[, 2]), c(data[, 3]))
  # create vector of x - axis coordinates
  labs <- c(1:length(data[, 1]))
  # create palindromic vector of y labels
  ylabs <- c(rev(interval_labels), interval_labels)
  limit <- length(interval_labels)
  breaks_limit <- limit - 0.5
  plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = labs)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        x = labs,
        xmin = as.numeric(labs) - 0.45,
        xmax = as.numeric(labs) + 0.45,
        ymin = c(data[, 2]),
        ymax = c(data[, 3])
      ),
      fill = "yellow2",
      alpha = 0.8
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        x = labs,
        xmin = as.numeric(labs) - 0.45,
        xmax = as.numeric(labs) + 0.45,
        ymin = c(data[, 4]) * -1,
        ymax = c(data[, 5]) * -1
      ),
      fill = "skyblue2",
      alpha = 0.8
    ) +
    ggplot2::scale_x_continuous(breaks = 1:length(labels), labels = labels) +
    ggplot2::scale_y_continuous(
      breaks = seq(-breaks_limit, breaks_limit, 1),
      limits = c(-limit, limit),
      labels = ylabs,
      minor_breaks = seq(-limit, limit, 0.5)
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = x_text_size,
        angle = 45,
        hjust = 1
      )
    ) +
    ggplot2::geom_hline(yintercept = 0,
               linetype = "dashed",
               size = 0.2) +
    ggplot2::geom_text(y = 0, ggplot2::aes(label = labels)) +
    ggplot2::coord_flip()
  return(plot_obj)
}
