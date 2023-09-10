\name{pois_PRC}
\alias{pois_PRC}

\title{PRC for Poisson data with rate parameter unknown
}

\description{
pois_PRC is used to derive the Predictive Ratio CUSUM (PRC) for Poisson data, when the rate parameter is unknown.
}

\usage{
pois_PRC( data = NULL, s = NULL, historical_data = NULL, historical_s = NULL,
          c0 = 1/2, d0 = 0, alpha_0 = NULL, k=2, two.sided=FALSE,
          h = log(100), FIR = FALSE, fFIR = 1/2, dFIR = 3/4,
          summary_list = TRUE, PRC_PLOT = TRUE, pdf_report = FALSE,
          path_pdf_report = tempdir(),
          xlab = "Observation Order", ylab = "PRC cumulative statistics",
          main = "PRC Poisson with unknown rate" )
}

\arguments{
  \item{data}{ vector; a univariate dataset for PCC implementation. Data needs to be in a vector form.
}
  \item{s}{ vector; refers to the period of time or grouping of counts, if data are rate dataset. It needs to be in a vector form.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{historical_s}{vector; refers to the period of time or grouping of counts, if historical_data are rate dataset. It needs to be in a vector form.
}
  \item{c0}{ scalar (positive); scalar (positive); hyperparameter of the prior G(c0, d0). It needs to be a number. The default is 1/2 and it refers to the initial reference prior G(1/2, 0).
}
  \item{d0}{ scalar (positive); scalar (positive); hyperparameter of the prior G(c0, d0). It needs to be a number. The default is 0 and it refers to the initial reference prior G(1/2, 0).
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{k}{ scalar (positive); tuning parameter, which represents deviation between IC and OOC state. Precisely, if k>1 then the OOC state corresponds to a (k-1)100\% increase in the rate parameter and if  k<1 corresponds (1-k)100\% decrease. The dafault value is 2.
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

pois_PRC provides the Predictive Ratio CUSUM (PRC) process for Poisson data, when the rate parameter parameter is unknown. The PRC process is based on sequential comparison between the IC and OOC state, via their predictive distributions. The parameter k represents the deviance between these two states. The default value is k=2, which corresponds to a doubling of the rate parameter. The testing starts as early as the second observation.

The initial prior for the unknown probability is a Gamma: G(c0, d0). Furthermore, the direct use of available historical data is possible via the power prior, if they are available. In this case, the default value for the power prior parameter alpha_0 is the reciprocal of the length of the historical data, which conveys the weight of a single observation to the prior information. The default prior is the non-informative reference prior G(1/2, 0), without the use of historical data.

A FIR option is available by inlfating the first few PRC statistics. The choice of the decision limit h reflects the false alarm tolerance.

}



\examples{
# the PRC process for the second application in
# "Design and properties of the Predictive Ratio Cusum (PRC) control charts"

### CURRENT DATA
CD <- c(1, 0, 0, 0, 1, 0, 3, 3, 3, 2, 5, 5, 2, 4, 4, 3, 4, 3, 8, 3, 2, 2)

### product exposures per million

sn <- c( 0.206, 0.313, 0.368, 0.678, 0.974, 0.927, 0.814, 0.696, 0.659, 0.775, 0.731,
         0.710, 0.705, 0.754, 0.682, 0.686, 0.763, 0.833, 0.738, 0.741, 0.843, 0.792 )

# regular process
pois_PRC(data = CD, s = sn)

# FIR process
pois_PRC(data = CD, s = sn, FIR = TRUE)

}
