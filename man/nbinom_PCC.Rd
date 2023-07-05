\name{nbinom_PCC}
\alias{nbinom_PCC}

\title{PCC for Negative Binomial data with probability parameter unknown
}

\description{
nbinom_PCC is used to derive the Predictive Control Chart (PCC) for Negative Binomial data, when the probability parameter is unknown and the number of failures is known.
}

\usage{
nbinom_PCC( data = NULL, historical_data = NULL,
            rl = NULL, a0 = 0, b0 = 1/2, alpha_0 = NULL,
            ARL_0 = 370.4, FAP = NULL, FIR = FALSE, fFIR = .99, aFIR = 1/8,
            summary_list = TRUE, PCC_PLOT = TRUE, PriorPosterior_PLOT = FALSE,
            historical_data_PLOT = FALSE, pdf_report = FALSE, path_pdf_report = tempdir(),
            xlab = "Observations", ylab = "Quality characteristic Values",
            main = "PCC Negative Binomial with unknown probability" )
}

\arguments{
  \item{data}{ vector; a univariate dataset for PCC implementation. Data needs to be in a vector form.
}
  \item{rl}{ scalar (positive); parameter of the likelihood. It is the known number of failures of the data and needs to be a number.
}
  \item{historical_data}{vector; an optional univariate dataset of historical data. Historical data needs to be in a vector form.
}
  \item{a0}{ scalar (positive); scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is 0 and it refers to the initial reference prior Beta(0, 1/2).
}
  \item{b0}{ scalar (positive); scalar (positive); shape hyperparameter of the prior Beta(a, b). It needs to be a number. The default is 1/2 and it refers to the initial reference prior Beta(0, 1/2).
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controlling the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{ARL_0}{scalar (positive); In Control (IC) Average Run Length (ARL). It is average number of IC data points that we will plot in the PCC before a false alarm occurs. The default value is 370.4
}
  \item{FAP}{scalar (between 0 and 1); False Alarm Probability (FAP). It is the probability of raising at least one false alarm out of a pre-determined number of N hypothesis tests and it can be used instead of ARL_0. It is based on the Sidak's correction.
}
  \item{FIR}{logical; If TRUE, then the Fast Initial Response (FIR) PCC is applied, which is an adjustment (narrows the PCC limits) for the initial fisrt tests.
}
  \item{fFIR}{a number between 0 and 1; It is used if FIR=TRUE. The default value is 0.99 and represents the proportion of the adjusted PCC region over the initial one for the first test.
}
  \item{aFIR}{non-negative number; It is used if FIR=TRUE. The default value is 0.125 and it is a smoothing parameter for the FIR adjustment.
}
  \item{summary_list}{logical; If it is TRUE, then a data frame is provided, containing the data sequence, the PCC limits along with the possible occurence of an alarm. It is TRUE by default.
}
  \item{PCC_PLOT}{logical; if TRUE, the PCC plot is displayed. It is TRUE by default.
}
  \item{PriorPosterior_PLOT}{ logical; if TRUE, the Prior and Posterior distribution(s) of the parameter(s) are displayed. It is TRUE by default.
}
  \item{historical_data_PLOT}{logical; if TRUE, the historical data are plotted precending the PCC plot. It is FALSE by default.
}
  \item{pdf_report}{logical; if TRUE then the summary list, PCC plot and Prior/Posterior distribution(s) of the parameter(s) are reported in a pdf file.
}
  \item{path_pdf_report}{specified path to save pdf report if pdf_report=TRUE. The default path is the temporary direction tempdir().
}
  \item{xlab, ylab, main}{The titles of the x-axis, y-axis and the overall title for the PCC plot. The default values are "Observation Order", "Quality characteristic Values" and "PCC Negative Binomial with unknown probability" respectively
}
}
\details{

nbinom_PCC provides the Predictive Control Chart (PCC) process for Negative Binomial data, when the probability parameter is unknown and the number of failures is known. The PCC process is based on the sequential testing of the future observable against the Highest Predictive Density (HPrD), which is obtained by the posterior predictive distribution. The PCC testing starts as early as the second observation.

The initial prior for the unknown probability is a Beta: Beta(a0, b0). Furthermore, the direct use of available historical data is possible via the power prior, if they are available. In this case, the default value for the power prior parameter alpha_0 is the reciprocal of the length of the historical data, which conveys the weight of a single observation to the prior information. The default prior is the non-informative reference prior Beta(0, 1/2), without the use of historical data.

A FIR option is available by narrowing the first few control limits. The metrics that can be used to control the false alarms tolerance are either ARL_0 or FAP.

}
\value{
Runs the PCC process, while if FIR=T, then the FIR-PCC is implemented. If plot=T, then a data frame is provided with the data, the upper and the lower limits of the Highest Predictive Density or Mass (HPrD/M) region, and an indication in case of alarms. If PCC_PLOT=T, then the PCC process is visualized, while if historical_data_PLOT=T, the historical data are added to the plot. If PriorPosterior_PLOT = T, a plot with the prior and the posterior is provided. Finally, if pdf_report=T, all the output is reported and saved into a pdf file.
}
\examples{
# 30 Negative Binomial observations introducing an outlier at the 15th observation
set.seed(5)
SimData <- rnbinom(n = 30, size = 10, prob = 0.9)
SimData[15] <- round( SimData[15] + 3*sqrt(10*0.1/(0.9^2)) )
nbinom_PCC(SimData, rl = 10)
}
