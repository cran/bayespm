\name{norm_HD}
\alias{norm_HD}

\title{The Highest Density (HD) interval of Normal distribution.
}
\description{
norm_HD is used to derive the narrowest interval of a predetermined coverage from a Normal distribution, containing the Highest Density (HD) values.
}

\usage{
norm_HD( cover = NULL, mu = NULL, sdv = NULL, plot = FALSE,
         xlab = "x", ylab = "Density" )
}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Density (HD). }
    \item{mu}{ scalar (real); It is the mean of the Normal distribution and needs to be a number.
}
  \item{sdv}{ scalar (positive); It is the standand deviation of the Normal distribution and needs to be a number.
}
  \item{plot}{logical; if TRUE, the plot of the HD region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HD plot. The default values are "x" and "Density".
}
}

\details{
norm_HD provides the Highest Density (HD) of a Normal distribution. HD is the narrowest interval with a pretermined coverage. The values included in HD have density larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
norm_HD(0.95, mu = 10, sdv = 1/2, plot = TRUE)

}
