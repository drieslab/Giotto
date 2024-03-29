% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience.R
\name{.createGiottoCosMxObject_all}
\alias{.createGiottoCosMxObject_all}
\title{Load and create a CosMx Giotto object from subcellular and aggregate info}
\usage{
.createGiottoCosMxObject_all(
  dir_items,
  FOVs,
  remove_background_polygon = TRUE,
  background_algo = c("range"),
  remove_unvalid_polygons = TRUE,
  cores,
  verbose = TRUE,
  instructions = NULL,
  ...
)
}
\arguments{
\item{dir_items}{list of full directory paths from \code{.read_cosmx_folder}}

\item{FOVs}{field of views to load (only affects subcellular data and images)}

\item{remove_background_polygon}{try to remove background polygon (default: FALSE)}

\item{background_algo}{algorithm to remove background polygon}

\item{remove_unvalid_polygons}{remove unvalid polygons (default: TRUE)}

\item{cores}{how many cores or threads to use to read data if paths are provided}

\item{verbose}{be verbose when building Giotto object}

\item{instructions}{list of instructions or output result from \code{\link[GiottoClass]{createGiottoInstructions}}}
}
\description{
Load and create a CosMx Giotto object from subcellular and aggregate info
}
\details{
Both \emph{subcellular} (subellular transcript detection information) and
\emph{aggregate} (aggregated detection count matrices by cell polygon from NanoString)
data will be loaded in. The two will be separated into 'cell' and 'cell_agg'
spatial units in order to denote the difference in origin of the two.
}
\seealso{
createGiottoCosMxObject .createGiottoCosMxObject_aggregate
.createGiottoCosMxObject_subcellular
}
\keyword{internal}
