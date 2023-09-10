\name{norm_mean2_PRC_h}
\alias{norm_mean2_PRC_h}

\title{Derivation of the decision limit for the PRC for Normal data with unknown parameters (mean)
}

\description{
norm_mean2_PRC_h is used to estimate the decision limit h for the Predictive Ratio CUSUM (PRC) for the mean of Normal data, when both the mean and the variance are unknown. The process controls the false alarms tolerance, providing the ARL_0 or the FAP, based on the standardized predictive distribution.
}
\usage{
norm_mean2_PRC_h( ARL_0 = 370.4, FAP = NULL, N=NULL, historical_data = NULL,
                  l0 = 0, a0 = -1/2, alpha_0 = NULL, k = 1, it = 1e4,
                  ARL0tol = 10/it )
}

\arguments{
   \item{ARL_0}{scalar (positive); lower bound of the In Control (IC) Average Run Length (ARL). It is average number of IC data points that we will plot in the PCC before a false alarm occurs. The default value is 370.4.
}
  \item{FAP}{scalar (between 0 and 1); lower bound of the False Alarm Probability (FAP). It is the probability of raising at least one false alarm out of a pre-determined number of N hypothesis tests and it can be used instead of ARL_0.
}
  \item{N}{ scalar (natural); number of N hypothesis tests. It is used for the FAP determination.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{l0}{ scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is 0.
}
  \item{a0}{ scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is -1/2.
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{k}{ scalar (positive); tuning parameter, which represents deviation between IC and OOC state. Precisely, if k>1 then the OOC state corresponds to a k standard deviations increase in the mean and if  k<1 corresponds to a k standard deviations decrease. The default value is 1.
}
  \item{it}{scalar (natural); iterations of the simulation. The default value is 10000.
}
  \item{ARL0tol}{scalar (positive); relative tolerance of ARL_0 estimate. It represents the maximum amount of error which is allowed. The default value is 10/it.
}
}
\details{



norm_mean2_PRC_h estimates the desicion limit h regarding the desired upper bound of false alarms tolerance for the Predictive Ratio CUSUM (PRC) process for the mean of Normal data, when both the mean and variance are unknown. The estimate concerns an one-sided PRC process and the default setup corresponds to one standard deviation increase in the mean. The estimate is based on the use of the standardized predictive distribution and the incorporation of historical data is possible via the power prior.

The default metric of the false alarms tolerance is ARL_0=370.4, while a FAP option is available. The desicion limit is obtained via the numerical method regula falsi when ARL_0 is selected, while simulated samples of fixed size are used when FAP is of interest.

}



\examples{

# Derivation of the decision limit of the first application in
# "Design and properties of the Predictive Ratio Cusum (PRC) control charts"

CD <- c( 31.0, 30.0, 32.0, 28.0, 33.2, 33.2, 35.1, 35.1, 33.9, 37.9,
         33.2, 36.5, 33.2, 35.1, 34.5, 36.5, 33.2, 35.1, 37.2, 32.6, 36.5 )
HD <- c( 31, 30, 33, 30, 33, 30, 31, 32, 32, 30, 33, 31, 34, 31, 34, 34, 36, 30,
         33, 29, 34, 32, 32, 28, 34, 32, 32, 30, 31, 29, 31, 29, 31, 32, 34, 34, 32 )

N <- length(CD)
n0 <- length(HD)
Pa0 <- 1/n0
M0F <- 31.8
LF <- 1/2
AF <- 2
BF <- 2.1^2
M0F ; LF ; AF ; BF

# To replicate results from application set 'it = 1e5'
norm_mean2_PRC_h( ARL_0 = NULL, FAP = 0.05, N = N, l0 = LF, a0 = AF,
                  historical_data = HD, alpha_0 = Pa0, it = 1e4 )





}
