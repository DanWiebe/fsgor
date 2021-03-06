% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_ancestry_annot.R
\name{find_boundaries}
\alias{find_boundaries}
\title{DAG roots and leaves extractor}
\usage{
find_boundaries(rellist)
}
\arguments{
\item{rellist}{list contains child terms as keys and nearest parents as values}
}
\value{
list contains root and leaf terms
}
\description{
DAG roots and leaves extractor
}
\examples{
\dontrun{
    onto_data <- fsgor::read_obo("~/go.obo")
    find_boundaries(onto_data$relation_data)
}
}
