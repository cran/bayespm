\name{t_HD}
\alias{t_HD}

\title{The Highest Density (HD) interval of Student's t distribution.
}
\description{
t_HD is used to derive the narrowest interval of a predetermined coverage from a Student's t distribution, containing the Highest Density (HD) values.
}
\usage{
t_HD( cover = NULL, df = NULL, mu = NULL, sdv = NULL, plot = FALSE,
      xlab = "x", ylab = "Density" )}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Density (HD). }
    \item{df}{ scalar (positive); It is the degrees of freedom of the Student's t distribution and needs to be a number.
}
    \item{mu}{ scalar (real); It is the mean of the Student's t distribution and needs to be a number. The default value is 0.
}
  \item{sdv}{ scalar (positive); It is the scale parameter of the Student's t distribution and needs to be a number. The default value is 1.
}
  \item{plot}{logical; if TRUE, the plot of the HD region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HD plot. The default values are "x" and "Density".
}
}

\details{
t_HD provides the Highest Density (HD) of a Student's t distribution. HD is the narrowest interval with a pretermined coverage. The values included in HD have density larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
t_HD( 0.95, df = 2, mu = 2, sdv = 3, plot = TRUE )
}
