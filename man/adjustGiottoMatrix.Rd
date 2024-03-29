% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary_giotto.R
\name{adjustGiottoMatrix}
\alias{adjustGiottoMatrix}
\title{Adjust expression values}
\usage{
adjustGiottoMatrix(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  batch_columns = NULL,
  covariate_columns = NULL,
  return_gobject = TRUE,
  update_slot = c("custom")
)
}
\arguments{
\item{gobject}{giotto object}

\item{spat_unit}{spatial unit (e.g. "cell")}

\item{feat_type}{feature type (e.g. "rna", "dna", "protein")}

\item{expression_values}{expression values to use}

\item{batch_columns}{metadata columns that represent different batch (max = 2)}

\item{covariate_columns}{metadata columns that represent covariates to regress out}

\item{return_gobject}{boolean: return giotto object (default = TRUE)}

\item{update_slot}{expression slot that will be updated (default = custom)}
}
\value{
giotto object or exprObj
}
\description{
Adjust expression values to account for known batch effects or technological covariates.
}
\details{
This function implements the \code{\link[limma]{removeBatchEffect}} function to
remove known batch effects and to adjust expression values according to provided covariates.
}
