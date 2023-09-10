\name{binom_PRC}
\alias{binom_PRC}

\title{PRC for Binomial data with probability parameter unknown
}

\description{
binom_PRC is used to derive the Predictive Ratio CUSUM (PRC) for Binomial data, when the probability parameter is unknown.
}

\usage{
binom_PRC( data = NULL, n = NULL, historical_data = NULL, historical_n = NULL,
           a0 = 1/2, b0 = 1/2, alpha_0 = NULL, k=2, two.sided=FALSE,
           h = log(100), FIR = FALSE, fFIR = 1/2, dFIR = 3/4,
           summary_list = TRUE, PRC_PLOT = TRUE, pdf_report = FALSE,
           path_pdf_report = tempdir(),
           xlab = "Observation Order", ylab = "PRC cumulative statistics",
           main = "PRC Binomial with unknown probability" )
}

\arguments{
  \item{data}{ vector; a univariate dataset for PCC implementation. Data needs to be in a vector form.
}
  \item{n}{ vector; number of trials for given dataset. It needs to be in a vector form.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{historical_n}{ vector; number of trials for given historical dataset. It needs to be in a vector form.
}
  \item{a0}{ scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is 1/2 and it refers to the initial reference prior Beta(1/2, 1/2).
}
  \item{b0}{ scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is 1/2 and it refers to the initial reference prior Beta(1/2, 1/2).
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{k}{ scalar (positive); tuning parameter, which represents deviation between IC and OOC state. Precisely, if k>1 then the OOC state corresponds to a (k-1)100\% increase in the expected odds of the success probability p and if  k<1 corresponds (1-k)100\% decrease. The dafault value is 2.
}
  \item{two.sided}{logical; If TRUE, then a two sided PRC is employed with tuning parameters k and 1/k respectively. It is FALSE by default.
}
  \item{h}{scalar (positive); Decision limit and we raise an alarm if the cumulative statistic exceeds it (or ?h for two sided shifts). The default value is log(100)=4.605, fow which we have a decisive cumulative evidence in favor of the OOC state.
}
  \item{FIR}{logical; If TRUE, then the Fast Initial Response (FIR) PRC is applied, which is an exponentially decreasing adjustment for the initial fisrt tests, inflating the PRC statistics.
}
  \item{fFIR}{non-negative number; It is used if FIR=TRUE. It represents the proportion of the inflation for the PRC statistic in the for the fisrt test. The default value is 1/2.
}
  \item{dFIR}{a number between 0 and 1; It is used if FIR=TRUE. It is a smoothing parameter, specifying the exponential decay of the adjustment (the smaller the d the fastest the decay). The default value is 3/4.
}
  \item{summary_list}{logical; If it is TRUE, then a data frame is provided, containing the data sequence, the PCC limits along with the possible occurence of an alarm. It is TRUE by default.
}
  \item{PRC_PLOT}{logical; if TRUE, the PRC plot is displayed. It is TRUE by default.
}
  \item{pdf_report}{logical; if TRUE then the summary list, PCC plot and Prior/Posterior distribution(s) of the parameter(s) are reported in a pdf file.
}
  \item{path_pdf_report}{specified path to save pdf report if pdf_report=TRUE. The default path is the temporary direction tempdir().
}
  \item{xlab, ylab, main}{The titles of the x-axis, y-axis and the overall title for the PCC plot. The default values are "Observation Order", "PRC cumulative statistics" and "PRC Binomial with unknown probability" respectively
}
}
\details{

binom_PRC provides the Predictive Ratio CUSUM (PRC) process for Binomial data, when the probability parameter is unknown. The PRC process is based on sequential comparison between the IC and OOC state, via their predictive distributions. The parameter k represents the deviance between these two states. The default value is k=2, which corresponds to a doubling of the expected odds ratio of the probability of success. The testing starts as early as the second observation.

The initial prior for the unknown probability is a Beta: Beta(a, b). Furthermore, the direct use of available historical data is possible via the power prior, if they are available. In this case, the default value for the power prior parameter alpha_0 is the reciprocal of the length of the historical data, which conveys the weight of a single observation to the prior information. The default prior is the non-informative reference prior Beta(1/2, 1/2), without the use of historical data.

A FIR option is available by inlfating the first few PRC statistics. The choice of the decision limit h reflects the false alarm tolerance.

}



\examples{

# the PRC process for the first 30 data points in the third application in
# "Design and properties of the Predictive Ratio Cusum (PRC) control charts"


### HISTORICAL DATA (FIRST BATCH)
HD <- c( 3, 3, 1, 5, 2, 1, 3, 3, 3, 0, 2, 1, 2, 1, 4, 1, 1, 0, 3, 2, 4, 6, 0, 1, 3, 2, 2, 4, 2, 1 )

### Bernoulli trials
nn <- 50

N0 <- length(HD)
NN0 <- rep(50, N0)

binom_PRC(data = HD, n = NN0)


# the PRC process for the last 10 data points in the third application in
# "Design and properties of the Predictive Ratio Cusum (PRC) control charts"

### prior parameters before the first batch

a0 <- 1/2
b0 <- 1/2

### posterior parameters after the first batch

ap <- sum(HD) + a0
bp <- sum(NN0) - sum(HD) + b0

hl = 4.332 # the decision limit is derived by the function binom_PRC_h

### CURRENT DATA (SECOND BATCH)
CD <- c(2, 4, 5, 2, 4, 8, 4, 4, 8, 5)

N <- length(CD)
NN <- rep(50, N)

binom_PRC(data = CD, n = NN, a0 = ap, b0 = bp, h = hl)

}
