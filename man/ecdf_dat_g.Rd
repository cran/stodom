% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stodom_functions.R
\name{ecdf_dat_g}
\alias{ecdf_dat_g}
\title{values of two ecdf and their cumulative difference}
\usage{
ecdf_dat_g(data_1, data_2, bins_size)
}
\arguments{
\item{data_1}{data 1.}

\item{data_2}{data 2.}

\item{bins_size}{bin size.}
}
\value{
The function returns a data table.
}
\description{
This function computes the values of two empirical cumulative distribution function as well as their cumulative differences.
}
\details{
This function computes the values of two empirical cumulative distribution function as well as their cumulative differences.
}
\examples{

# load stodom
require(stodom)

 data_a <- rnorm(500, 3, 2)
 data_b <- rnorm(500, 1, 2)

# compute the values of two ecdfs and their cumulative differences.
ecdf_dat_g(data_1 = data_a, data_2 = data_b, bins_size = 1)
}
