\name{pois_PRC_h}
\alias{pois_PRC_h}

\title{Derivation of the decision limit for the PRC for Poisson data with probability parameter unknown
}

\description{
pois_PRC_h is used to estimate the decision limit h for the Predictive Ratio CUSUM (PRC) for Poisson data. The process controls the false alarms tolerance, providing a lower bound of ARL_0 or FAP, based on the  marginal distribution.
}
\usage{
pois_PRC_h( ARL_0 = 370.4, FAP = NULL, N = NULL, s = NULL,
            historical_data = NULL, historical_s = NULL,
            c0 = NULL, d0 = NULL, alpha_0 = NULL, k = 2,
            it = 1e4, ARL0tol = 10/it )
}

\arguments{
   \item{ARL_0}{scalar (positive); lower bound of the In Control (IC) Average Run Length (ARL). It is average number of IC data points that we will plot in the PCC before a false alarm occurs. The default value is 370.4.
}
  \item{FAP}{scalar (between 0 and 1); lower bound of the False Alarm Probability (FAP). It is the probability of raising at least one false alarm out of a pre-determined number of N hypothesis tests and it can be used instead of ARL_0.
}
  \item{N}{ scalar (natural); number of N hypothesis tests. It is used for the FAP determination.
}
  \item{s}{ vector; refers to the period of time or grouping of counts, if data are rate dataset. It needs to be in a vector form.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{historical_s}{vector; refers to the period of time or grouping of counts, if historical_data are rate dataset. It needs to be in a vector form.
}
  \item{c0}{ scalar (positive); scalar (positive); hyperparameter of the prior G(c0, d0). It needs to be a number.
}
  \item{d0}{ scalar (positive); scalar (positive); hyperparameter of the prior G(c0, d0). It needs to be a number.
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{k}{ scalar (positive); tuning parameter, which represents deviation between
IC and OOC state. Precisely, if k>1 then the OOC state corresponds to a (k-1)100\% increase in the rate parameter and if  k<1 corresponds (1-k)100\% decrease. The dafault value is 2.
}
  \item{it}{scalar (natural); iterations of the simulation. The default value is 10000.
}
  \item{ARL0tol}{scalar (positive); relative tolerance of ARL_0 estimate. It represents the maximum amount of error which is allowed. The default value is 10/it.
}
}
\details{



pois_PRC_h estimates the desicion limit h regarding the desired upper bound of false alarms tolerance for the Predictive Ratio CUSUM (PRC) process for Binomial data, when the probability parameter is unknown. The estimate concerns an one-sided PRC process and the default setup corresponds to a doubling of the rate parameter (k=2). The estimate is based on the use of the  marginal distribution (prior predictive) and the incorporation of historical data is possible via the power prior.

The default metric of the false alarms tolerance is ARL_0=370.4, while a FAP option is available. The desicion limit is obtained via the numerical method regula falsi when ARL_0 is selected, while simulated samples of fixed size are used when FAP is of interest.

}



\examples{



pois_PRC_h(ARL_0 = 150, c0 = 40, d0 = 10, it = 1e3)


}
