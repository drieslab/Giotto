% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience.R
\name{.visium_read_scalefactors}
\alias{.visium_read_scalefactors}
\title{Read Visium ScaleFactors}
\usage{
.visium_read_scalefactors(json_path = NULL)
}
\arguments{
\item{json_path}{path to scalefactors_json.json for Visium experimental data}
}
\value{
scalefactors within the provided json file as a named list,
or NULL if not discovered
}
\description{
Read Visium ScaleFactors
}
\details{
asserts the existence of and reads in a .json file
containing scalefactors for Visium data in the expected format.
Returns NULL if no path is provided or if the file does not exist.
}
\keyword{internal}
