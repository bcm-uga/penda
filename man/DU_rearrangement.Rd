% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/choose_param.R
\name{DU_rearrangement}
\alias{DU_rearrangement}
\title{Reorganization of results matrix after the multiple test.}
\usage{
DU_rearrangement(multiple_tests, simu_data, multiple_values)
}
\arguments{
\item{multiple_tests}{The results matrix, which has a form of "DUDUDUDU".}

\item{simu_data}{The matrix of initial_data or simulated_data in the simulation, used for dimnames.}

\item{multiple_values}{The vector with different values of the test.}
}
\value{
This function returns a list of two logical matrices. The D matrix with TRUE
for genes down-regulated, and the U matrix with TRUE for genes up-regulated.
}
\description{
This function reorganizes matrices resulting from the multiple test. It takes a
"DUDUDUDU" matrix and make two independant D and U matrices.
}
