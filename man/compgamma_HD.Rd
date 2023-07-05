\name{compgamma_HD}
\alias{compgamma_HD}

\title{The Highest Density (HD) interval of Compound Gamma distribution.
}
\description{
compgamma_HD is used to derive the narrowest interval of a predetermined coverage from a Compound Gamma distribution, containing the Highest Density (HD) values.
}
\usage{
compgamma_HD( cover = NULL, a = NULL, c = NULL, d = NULL, plot = FALSE,
              xlab = "x", ylab = "Density" )}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Density (HD). }
    \item{a}{ scalar (positive); It is the shape1 parameter of the Compound Gamma distribution and needs to be a number.
}
    \item{c}{ scalar (positive); It is the shape2 parameter of the Compound Gamma distribution and needs to be a number.
}
  \item{d}{ scalar (positive); It is the scale parameter of the Compound Gamma distribution and needs to be a number.
}
  \item{plot}{logical; if TRUE, the plot of the HD region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HD plot. The default values are "x" and "Density".
}
}

\details{
compgamma_HD provides the Highest Density (HD) of a Compound Gamma distribution. HD is the narrowest interval with a pretermined coverage. The values included in HD have density larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
compgamma_HD(0.95, 2, 10, 10, plot = TRUE)

}
