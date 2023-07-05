\name{nbinom_HM}
\alias{nbinom_HM}

\title{The Highest Mass (HM) interval of Beta-Negative Binomial distribution.
}
\description{
nbinom_HM is used to derive the narrowest interval of a predetermined coverage from a Negative Binomial distribution, containing the Highest Mass (HM) values.
}
\usage{
nbinom_HM( cover = NULL, r = NULL, p = NULL, plot = FALSE,
           xlab = "x", ylab = "Probability" )}
\arguments{
  \item{cover}{ scalar (between 0 and 1); The coverage (probability) of the Highest Mass (HM).}
    \item{r}{ scalar (integer); It is the number of failures until the experiment is stopped.
}
  \item{p}{ scalar (between 0 and 1); It is the probability of success for each Bernoulli trial.
}
  \item{plot}{logical; if TRUE, the plot of the HM region is displayed. It is FALSE by default.
}
  \item{xlab, ylab}{The titles of the x-axis and y-axis for the HM plot. The default values are "x" and "Probability".
}
}

\details{
nbinom_HM provides the Highest Mass (HM) of a Negative Binomial distribution. HM is the narrowest region, with the minimum absolute difference from a pretermined coverage. The values included in HM have probability larger than every point outside the region.
}
\value{
Returns a data frame with the upper and the lower limits of the Highest Predictive Denisty or Mass (HPrD/M) region, the lower and the achieved coverage. If plot=T, then the HPrD/M is visualized.
}
\examples{
nbinom_HM(0.95, 4, 0.2, plot = TRUE)

}
