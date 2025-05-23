% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTest.R
\name{CTest}
\alias{CTest}
\title{CT estimation of expct}
\usage{
CTest(
  differentialtimevaryingpredictors = differentialtimevaryingpredictors,
  outcome = outcome,
  namesofnewpredictorvariables = namesofnewpredictorvariables,
  laglongreducedummy = laglongreducedummy,
  method = "bam",
  gamma = gamma,
  numberofknots = numberofknots,
  ktrend = 3,
  debug = debug,
  lengthcovariates = lengthcovariates,
  plot_show = FALSE,
  weights = NULL
)
}
\arguments{
\item{differentialtimevaryingpredictors}{The variables that will be a varying-coefficient of differential time (AKA the lags you want to know what times they predict the outcome). This must be specified as a vector using c("variables here"). e.g. c("X","Y") (REQUIRED)}

\item{outcome}{This is each of the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2") (REQUIRED)}

\item{namesofnewpredictorvariables}{This is the name of the predictors.}

\item{laglongreducedummy}{This is the long data output from the data manipulation.}

\item{method}{Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".}

\item{gamma}{This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)}

\item{numberofknots}{The number of k selection points used in the model for stage 1 (see ?choose.k in mgcv package for more details) (note that this is for the raw data k2 refers to the k for the re-blocked data), default is 10. The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL)}

\item{ktrend}{The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)}

\item{debug}{This will print more useless information as it goes along. Only useful for troubleshooting problems. (OPTIONAL, UNCOMMONLY SPECIFIED)}

\item{lengthcovariates}{the number of covariates (+2)}

\item{plot_show}{Controls if show the plot}

\item{weights}{optional vector of weights to be passed to gam/bam}

\item{predictionstart}{The differential time value to start with, default is NULL, and the lowest time difference in the time series will be used (use lower value if you're first value if you're interested in a smaller interval prediction) e.g. predictionstart = 1. If this is not specified and using a continuous time model, make sure to set blockdata = TRUE so that it will be automatically chosen. (OPTIONAL)}

\item{predictionsend}{The differential time value to end with. This means how long you want your largest time difference in the study to be (i.e. if you wanted to predict up to allow time predictions up to 24 hours and your time intervals were specified in hours, you would set predictionsend = 24). If this is not specified and using a continuous time model, make sure to set blockdata = TRUE so that it will be automatically chosen. (OPTIONAL)}
}
\value{
The output of this function is: The estimated model
}
\description{
PLEASE USE THE expct FUNCTION RATHER THAN THIS UNLESS YOU WOULD LIKE TO SPECIFY THIS MANUALLY.
}
