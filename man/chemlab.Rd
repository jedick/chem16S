\encoding{UTF-8}
\name{chemlab}
\alias{chemlab}
\title{Formatted labels for chemical metrics}
\description{
  Formatted labels for indicating chemical metrics on plots.
}

\usage{
  chemlab(varname)
}

\arguments{
  \item{varname}{character, the name of a variable}
}

\details{
This function provides formatted labels for chemical metrics and related variables.
\code{varname} can be any one of the chemical metrics described in \code{\link{calc_metrics}}.
}

\seealso{
This function is used by \code{\link{plot_metrics}}, \code{\link{plot_ps_metrics}}, and \code{\link{plot_ps_metrics2}}.
}

\examples{
chemlab("Zc")
}