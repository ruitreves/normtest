% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odds.R
\name{combine_fc_pvals}
\alias{combine_fc_pvals}
\title{combine fold changes and pvals into single dataframe. takes column
i from data1 and column i from data2 and puts them next to each other}
\usage{
combine_fc_pvals(data1, data2)
}
\arguments{
\item{data1}{a data object with n columns and m rows}

\item{data2}{a data object with n columns and m rows}
}
\value{
a dataframe made of data1 and data2
}
\description{
combine fold changes and pvals into single dataframe. takes column
i from data1 and column i from data2 and puts them next to each other
}
