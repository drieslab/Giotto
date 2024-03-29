% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{create_jackstrawplot}
\alias{create_jackstrawplot}
\title{create_jackstrawplot}
\usage{
create_jackstrawplot(
  jackstraw_data,
  ncp = 20,
  ylim = c(0, 1),
  threshold = 0.01
)
}
\arguments{
\item{jackstraw_data}{result from jackstraw function (`testresult$p`)}

\item{ncp}{number of principal components to calculate}

\item{ylim}{y-axis limits on jackstraw plot}

\item{threshold}{p.value threshold to call a PC significant}
}
\value{
ggplot
}
\description{
create jackstrawplot with ggplot
}
\keyword{internal}
