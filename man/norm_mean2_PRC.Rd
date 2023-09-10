\name{norm_mean2_PRC}
\alias{norm_mean2_PRC}

\title{PRC for Normal data with unknown parameters (mean)
}

\description{
norm_mean2_PRC is used to derive the Predictive Ratio CUSUM (PRC) for the mean of Normal data, when both the mean and the variance are unknown.
}

\usage{
norm_mean2_PRC( data = NULL, historical_data = NULL, mu0 = 0, l0 = 0,
                a0 = -1/2, b0 = 0, alpha_0 = NULL, k = 1, two.sided = FALSE,
                h = log(100), FIR = FALSE, fFIR = 1/2, dFIR = 3/4,
                summary_list = TRUE, PRC_PLOT = TRUE, pdf_report = FALSE,
                path_pdf_report = tempdir(), xlab = "Observation Order",
                ylab = "PRC cumulative statistics",
                main = "PRC Normal with unknown parameters (mean model)" )
}

\arguments{
  \item{data}{ vector; a univariate dataset for PCC implementation. Data needs to be in a vector form.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{mu0}{ scalar; hyperparameter of the prior NIG(mu0, l0, a0, b0). It needs to be a number. The default is 0 and it refers to the initial reference prior NIG(0, 0, -1/2, 0).
}
  \item{l0}{ scalar (positive); hyperparameter of the prior NIG(mu0, l0, a0, b0). It needs to be a number. The default is 0 and it refers to the initial reference prior NIG(0, 0, -1/2, 0).
}
  \item{a0}{ scalar (positive); hyperparameter of the prior NIG(mu0, l0, a0, b0). It needs to be a number. The default is -1/2 and it refers to the initial reference prior NIG(0, 0, -1/2, 0).
}
  \item{b0}{ scalar (positive); hyperparameter of the prior NIG(mu0, l0, a0, b0). It needs to be a number. The default is 0 and it refers to the initial reference prior NIG(0, 0, -1/2, 0).
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{k}{ scalar (positive); tuning parameter, which represents deviation between IC and OOC state. Precisely, if k>1 then the OOC state corresponds to a k standard deviations increase in the mean and if  k<1 corresponds to a k standard deviations decrease. The default value is 1.
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

norm_mean2_PRC provides the Predictive Ratio CUSUM (PRC) process for the mean of Normal data, when both the mean and the variance are unknown. The PRC process is based on sequential comparison between the IC and OOC state, via their predictive distributions. The parameter k represents the deviance between these two states. The default value is k=1, which corresponds to one standard deviation increase in the mean. The testing starts as early as the second observation, apart from when using the default non-informative prior, without historical data. In this special case, the PCC testing starts at the third observation.

The initial prior for the unknown probability is a Normal Inverse Gamma: NIG(mu0, l0, a0, b0). Furthermore, the direct use of available historical data is possible via the power prior, if they are available. In this case, the default value for the power prior parameter alpha_0 is the reciprocal of the length of the historical data, which conveys the weight of a single observation to the prior information. The default prior is the non-informative reference prior NIG(0, 0, -1/2, 0), without the use of historical data.

A FIR option is available by inlfating the first few PRC statistics. The choice of the decision limit h reflects the false alarm tolerance.

}



\examples{


# the PRC process for the first application in
# "Design and properties of the Predictive Ratio Cusum (PRC) control charts"


### CD: Current data (New reagent)
### HD: Historical data (Previous reagent)

CD <- c( 31.0, 30.0, 32.0, 28.0, 33.2, 33.2, 35.1, 35.1, 33.9, 37.9,
         33.2, 36.5, 33.2, 35.1, 34.5, 36.5, 33.2, 35.1, 37.2, 32.6, 36.5 )
HD <- c( 31, 30, 33, 30, 33, 30, 31, 32, 32, 30, 33, 31, 34, 31, 34, 34, 36, 30,
         33, 29, 34, 32, 32, 28, 34, 32, 32, 30, 31, 29, 31, 29, 31, 32,34,34,32 )

N <- length(CD)
n0 <- length(HD)

### initial prior parameters

M0F <- 31.8
LF <- 1/2
AF <- 2
BF <- 2.1^2

norm_mean2_PRC( data = CD, historical_data = HD, alpha_0 = 1/n0, mu0 = M0F,
                l0 = LF, a0 = AF, b0 = BF, h = 3.749, two.sided = TRUE )



### a real data application to aPTT values

### CURRENT DATA aPTT
CD <- c( 29.0, 29.1, 28.7, 28.2, 28.0, 29.1, 28.6, 28.7, 28.6, 29.0, 28.4,
         28.1, 28.8, 29.7, 28.8, 29.8, 28.8, 29.4, 28.4, 28.7, 28.7, 29.5,
         28.5, 28.4, 28.1, 28.6, 28.2, 29.6, 28.9, 29.1, 29.0, 29.9, 28.6,
         29.3, 28.2, 28.6, 27.6, 27.3, 28.7, 27.2, 28.4, 28.0, 28.4, 27.8,
         28.4, 28.4, 27.7, 29.2, 27.5, 27.7)

### HISTORICAL DATA aPTT
HD <- c( 28.0, 28.9, 27.7, 29.3, 28.9, 29.5, 28.2, 27.5, 28.8, 28.9, 28.7,
         27.4, 28.6, 28.5, 29.6, 28.7, 21.3, 29.4, 28.1, 28.9, 28.3, 27.6,
         29.0, 29.2, 27.8, 29.1, 28.9, 29.4, 29.4, 28.9, 28.9, 29.2, 29.4,
         29.4, 28.1, 28.5, 29.7, 29.3, 28.6, 29.2, 29.3, 29.3, 29.3, 30.0,
         29.1, 29.1, 26.8, 29.0, 29.3, 28.3)


norm_mean2_PRC( data = CD, historical_data = HD, mu0 = 28.9,
                l0 = 1/4, a0 = 2, b0 = 0.49, two.sided = TRUE )


}
