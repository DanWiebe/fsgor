% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_data_viz_preparation.R
\name{create_intervals_matrix}
\alias{create_intervals_matrix}
\title{Create rectangle plot coordinates for each Go term}
\usage{
create_intervals_matrix(intervals)
}
\arguments{
\item{intervals}{- vector of intervals in character representation (e.g. "1-6")}
}
\value{
- matrix contains interval start coordinates in the first column and interval end coordinates
          in the second row (note: the start coordinate is deducted by 1 in order to make the appropriate layout for rectangles)
}
\description{
Create rectangle plot coordinates for each Go term
}
\examples{
\dontrun{
    create_intervals_matrix(intervals)
}
}
