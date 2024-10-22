% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova_tukey.R
\name{run_tukey}
\alias{run_tukey}
\title{Runs a Tukey multicomparison test on each row of my_data based on var1}
\usage{
run_tukey(my_data, var1)
}
\arguments{
\item{my_data}{a data object}

\item{var1}{a list that points column names of my_data to the factors of the experiment
return a data object containing p values for each groupwise comparison as specified by var1, computed by TukeyHSD}
}
\description{
Runs a Tukey multicomparison test on each row of my_data based on var1
}
