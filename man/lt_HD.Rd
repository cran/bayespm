\name{lt_HD}
\alias{lt_HD}

\title{The Highest Density (HD) interval of Logt distribution.
}
\description{
lt_HD is used to derive the narrowest interval of a predetermined coverage from a Logt distribution, containing the Highest Density (HD) values.
}
\usage{
lt_HD( cover = NULL, df = NULL, mulog = NULL, sdlog = NULL, plot = FALSE,
       xlab = "x", ylab = "Density" )}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Density (HD). The default value is 0.95.}
    \item{df}{ scalar (positive); It is the degrees of freedom of the Logt distribution and needs to be a number.
}
    \item{mulog}{ scalar (real); It is the scale parameter of the Logt distribution and needs to be a number. The default value is 0.
}
  \item{sdlog}{ scalar (positive); It is the shape parameter of the Logt distribution and needs to be a number. The default value is 1.
}
  \item{plot}{logical; if TRUE, the plot of the HD region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HD plot. The default values are "x" and "Density".
}
}

\details{
lt_HD provides the Highest Density (HD) of a Logt distribution. HD is the narrowest interval with a pretermined coverage. The values included in HD have density larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
lt_HD(0.95, 10, 0, 1/2, plot = TRUE)

}
