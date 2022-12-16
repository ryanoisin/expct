# to dos:
# - suppress detrending? k3 = FALSE? k3 = 0? check CTest (line 50) %%Solution: by adding an argument time_trend, if time_trend = FALSE, k3 = 0 %%FINISHED
# - make evertyhing into one function call %%Solution: add a boot argument in the CT_LAG function.
# - make version control on github  %%FINISHED
# - output looking the same %% FINISHED
# - efficiency (nick suggestions) %% Add an argument bam to indicate if we use the bam function to run these code.%% FINISHED
# - test bootstrap feasibility



#' Continuous Time-Varying Effect model (expct)
#'
#' This is the PRIMARY expct function which works for providing some lag information about continuous data.
#' @param dataset Specify the data frame that contains the interested variables, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
#' @param Time The name of the Time column in the data E.G. Time = "Time" (must be specified).
#' @param outcome This is the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2"). If it is NULL, it will consider each variables as outcome once.
#' @param ID The name of the ID column in the data E.G. ID = "ID"
#' @param estimate The relationship which we are interested, estimate = "marginal" or "partial". Default is the "marginal".
#' @param Tpred A vector which indicates that interested time points, e.g. seq(0,30,1)
#' @param plot_show The option to suppress the plot outcomes. The default if FALSE which means the plot outcomes will not appear.
#' @param boot Indicate if we perform bootstrapping estimation or not. If boot == True, we perfomr bootstrapping estimation. The default value is False
#' @param output_type Indicate which output form will be returned. If output_type == "CI", point estimations and corresponding CIs will be returned. If output_type =="PE", only ponit estimation will be returned. If output_type =="SCI", the Simultaneous CIs will be returned. The default value is "CI"
#' @param standardized This specifies whether all of the variables (aside from Time) should be standardized. Options are TRUE, FALSE, and "center". TRUE means within-person standardize each variable (aka get the person-centered z-scores), FALSE means use the raw data, "center" means to only within-person mean-center the variables. Default = TRUE. FALSE is not recommended unless you have done these transformations yourself (OPTIONAL)
#' @param method Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".
#' @param gamma This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param k The number of k selection points used in the model (see ?choose.k in mgcv package for more details) The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL, BUT RECOMMENDED)
#' @param ktrend The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)
#' @param time_trend This argument is used to suppress time trend or not. If time_trend = FALSE, k3 = 0. The default value is TRUE
#' @param iterations How many times bootstrapping single estimation you want to perform, the default value is 50
#' @param quantiles The quantiles to build bootstrapping CI, the default value is c(low_quantile, high_quantile) = c(.025, 0.975)
#' @param pivot Indicate the pivot of the bootstrapping CI is the mean of bootstrapping point estimations or the median of bootstrapping point estimations. The default value is ''Mean''
#' @param ncores How many cores you want to use. If it is null, ncores = detectCores()/2
#' @return The output of this function is: The point estimation of all specified marginal/partial effects (contained in a list).  If SE is true, all corresponding High-CIs and Low-CIs will also be returned (contained in a list).
#' @import mgcv
#' @import plyr
#' @import zoo
#' @import reshape2
#' @import Rcpp
#' @import foreach
#' @import parallel
#' @import pbapply
#' @import iterators
#' @export


expct <- function(dataset = NULL,
                  Time = "Time",
                  outcome = NULL,
                  ID = "ID",
                  estimate = "marginal",
                  Tpred = seq(0,30,1),
                  plot_show = FALSE ,
                  boot = FALSE,
                  output_type = "CI",
                  standardized = TRUE,
                  method = "bam",
                  gamma = 1,
                  k = 30,
                  ktrend = 3,
                  time_trend = TRUE,
                  iterations = 10,
                  quantiles = c(.025, 0.975),
                  pivot = "Mean",
                  #datamanipu = "DT",
                  ncores = NULL,
                  ctype = "PSOCK",
                  weighting = FALSE,
                  MBB_block = "non-Fixed",
                  ...

) {

  #LOAD NECESSARY PACKAGES
  #library(mgcv) #USED FOR THE PRIMARY ANALYSES
  #library(plyr)
  #library(zoo)
  #library(reshape2)
  #
  if(time_trend==FALSE){ # Suppress the time trend by setting k3 = 0
    ktrend = 0
  }

  if(boot == FALSE){
    Result = expct_single(
      dataset = dataset,
      Time = Time,
      outcome = outcome,
      ID = ID,
      estimate = estimate,
      Tpred = Tpred,
      plot_show = plot_show,
      quantiles = quantiles,
      boot = boot,
      output_type = output_type,
      standardized = standardized,
      method = method,
      gamma = gamma,
      k = k,
      #datamanipu = datamanipu,
      ktrend = ktrend,
      weighting = weighting,
      ...
    )
  }else if(boot == TRUE | boot  == "MBB" ){
    Result = expct_boot(
      dataset = dataset,
      Time = Time,
      outcome = outcome,
      ID = ID,
      estimate = estimate,
      iterations = iterations,
      quantiles = quantiles,
      boot = boot,
      pivot = pivot,
      Tpred = Tpred,
      plot_show = plot_show,
      ncores = ncores,
      standardized = standardized,
      method = method,
      gamma = gamma,
      k = k,
      ktrend = ktrend,
      ctype = ctype,
      MBB_block = MBB_block,
      ...
    )
  }

  return(Result)

}




