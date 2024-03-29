% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_interaction_spot.R
\name{do_permuttest_spot}
\alias{do_permuttest_spot}
\alias{.do_permuttest_original_spot}
\alias{.do_permuttest_random_spot}
\alias{.do_multi_permuttest_random_spot}
\alias{.do_permuttest_spot}
\title{Spot permutation testing}
\usage{
.do_permuttest_original_spot(
  sel_int,
  other_ints,
  select_ind,
  other_ind,
  name = "orig",
  proximityMat,
  expr_residual
)

.do_permuttest_random_spot(
  sel_int,
  other_ints,
  select_ind,
  other_ind,
  name = "perm_1",
  proximityMat,
  expr_residual,
  set_seed = TRUE,
  seed_number = 1234
)

.do_multi_permuttest_random_spot(
  sel_int,
  other_ints,
  select_ind,
  other_ind,
  proximityMat,
  expr_residual,
  n = 100,
  cores = NA,
  set_seed = TRUE,
  seed_number = 1234
)

.do_permuttest_spot(
  sel_int,
  other_ints,
  select_ind,
  other_ind,
  proximityMat,
  expr_residual,
  n_perm = 100,
  adjust_method = "fdr",
  cores = 2,
  set_seed = TRUE,
  seed_number = 1234
)
}
\description{
Test spot interactions using permutations
}
\section{Functions}{
\itemize{
\item \code{.do_permuttest_original_spot()}: Calculate original values for spots

\item \code{.do_permuttest_random_spot()}: Calculate random values for spots

\item \code{.do_multi_permuttest_random_spot()}: Calculate multiple random values for spots

\item \code{.do_permuttest_spot()}: Performs permutation test on subsets of a matrix for spots

}}
\keyword{internal}
