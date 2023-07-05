\name{betanbinom_HM}
\alias{betanbinom_HM}

\title{The Highest Mass (HM) interval of Beta-Negative Binomial distribution.
}
\description{
betanbinom_HM is used to derive the narrowest interval of a predetermined coverage from a Beta-Negative Binomial distribution, containing the Highest Mass (HM) values.
}
\usage{
betanbinom_HM( cover = NULL, r = NULL, a = NULL, b = NULL, plot = FALSE,
               xlab = "x", ylab = "Probability" )}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Mass (HM). }
    \item{r}{ scalar (integer); It is the number of failures of the Negative Binomial distribution and needs to be a number.
}
  \item{a}{ scalar (positive); It is the shape1 parameter of the Beta distribution and needs to be a number.
}
  \item{b}{ scalar (positive); It is the shape2 parameter of the Beta distribution and needs to be a number.
}
  \item{plot}{logical; if TRUE, the plot of the HM region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HM plot. The default values are "x" and "Probability".
}
}

\details{
betanbinom_HM provides the Highest Mass (HM) of a Beta-Negative Binomial distribution. HM is the narrowest region, with the minimum absolute difference from a pretermined coverage. The values included in HM have probability larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
betanbinom_HM(0.95, 5, 20, 80, plot = TRUE)

}
