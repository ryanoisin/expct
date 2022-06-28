
#' Continuous Time-Varying Effect model for a bootstrapping estimation.
#'
#' This is the CTVEM function which performs bootstrapping estimation. The output of this function is the point estimations, high and low cooresponding CIs.
#' @param data Specify the data frame that contains the interested data, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
#' @param Time The name of the Time variable. E.G. Time = "Time" (must be specified).
#' @param outcome This is the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2"). If it is NULL, it will consider each variables as outcome once.
#' @param ID The name of the ID column in the data E.G. ID = "ID"
#' @param estimate The relationship which we are interested, estimate = "marginal" or "partial". Default is the "marginal".
#' @param iterations How many times bootstrapping single estimation you want to perform, the default value is 50
#' @param quantiles The quantiles to build bootstrapping CI, the default value is c(low_quantile, high_quantile) = c(.025, 0.975)
#' @param boot Indicate if we perform bootstrapping estimation or not. If boot == True, we perfomr bootstrapping estimation. The default value is False
#' @param pivot Indicate the pivot of the bootstrapping CI is the mean of bootstrapping point estimations or the median of bootstrapping point estimations. The default value is ''Mean''
#' @param Tpred A vector which indicates that interested time points, e.g. seq(0,30,1)
#' @param plot_show The option to suppress the plot outcomes. The default if FALSE which means the plot outcomes will not appear.
#' @param ncores How many cores you want to use. If it is null, ncores = detectCores()/2
#' @param standardized This specifies whether all of the variables (aside from Time) should be standardized. Options are TRUE, FALSE, and "center". TRUE means within-person standardize each variable (aka get the person-centered z-scores), FALSE means use the raw data, "center" means to only within-person mean-center the variables. Default = TRUE. FALSE is not recommended unless you have done these transformations yourself (OPTIONAL)
#' @param method Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".
#' @param gamma This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param k The number of k selection points used in the model for stage 1 (see ?choose.k in mgcv package for more details) The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL, BUT RECOMMENDED)
#' @param ktrend The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)
#  return  a) point estimates, b) highCI, c) lowCI, but now based on the bootstrap

CTVEM_boot <-
  function(data = NULL,
           Time = "Time",
           outcome = NULL,
           ID = "ID",
           estimate = "marginal",
           iterations = 50,
           quantiles = c(.025, 0.975),
           boot = TRUE,
           pivot = "Mean",
           Tpred = seq(0, 30, 1),
           plot_show = FALSE,
           ncores = NULL,
           standardized = TRUE,
           method = "bam",
           gamma = 1,
           k = 10,
           ktrend = 3
           ) {
    if(is.null(ncores)){
    ncores = detectCores()/2
  }
  cl =  makeCluster(ncores)
  registerDoSNOW(cl)
  print(paste("Perform bootstrapping estimation with bootstrapping times = ", iterations, " ; Use ", ncores, " CPU cores"))
  pb = progress_bar$new(
    format = ":letter [:bar] :elapsed | eta: :eta",
    total = iterations,
    width = 60)
  progress_letter <- rep("Running", iterations)
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  }
  opts = list(progress = progress)
  bootstrap_results = foreach(iii=1:iterations, .combine="rbind",.options.snow = opts,.export=c("CTVEM_single","datamanipulation","CTest"),.packages = c("plyr","zoo","reshape2","mgcv","matrixStats")) %dopar%{

    Select = sort(sample(seq(1,nrow(data),1),size = nrow(data),replace = T),decreasing = F)
    data_select = data[Select,]
    data_select = data.frame(data_select)
    CTVEM_single(
      data = data_select,
      Time = Time,
      ID = ID,
      estimate = estimate,
      Tpred = Tpred,
      plot_show = plot_show,
      outcome = outcome,
      boot = boot,
      standardized = standardized,
      method = method,
      gamma = gamma,
      k = k,
      ktrend = ktrend
    )
  }
  stopCluster(cl)

  # Analyze bootstrapping results
  # Get how many time-varying effects we have
  colnames_data = colnames(data)
  varnames = colnames_data[-c(which(colnames_data == ID), which(colnames_data == Time))]
  if(is.null(outcome)){
    Result_length = length(varnames)*length(varnames)
  }else{ # with specified outcomes. E.g. three variables X1, X2 and X3. Outcome is c("X1","X2"). We only consider 3X2 = 6 marginal/partial effects
    Result_length = length(outcome)*length(varnames)
  }

  # Build lists to contain Bootstrapping point estimations and corresponding CIs

  Boot_preds = vector("list",length = Result_length) # Contain each marginal effects estimation
  Boot_highCI = vector("list",length = Result_length) # Contain Upper CI for each marginal effects estimation
  Boot_lowCI = vector("list", length = Result_length) # Contain Lower CI for each marginal effects estimation

  # Compute Bootstrapping point estimations and corresponding CIs

  for (nn in 1:Result_length){
    Bootstrap_Pre = bootstrap_results[seq(nn, nrow(bootstrap_results),Result_length),] # Extract all bootstrapping estimation one by one
    if(pivot == "Mean"){
      Boot_preds[[nn]] = colMeans(Bootstrap_Pre)
    }else{
      Boot_preds[[nn]] = colMedians(Bootstrap_Pre)
    }
    H_CI = c()
    L_CI = c()
    for (i in 1:ncol(Bootstrap_Pre)) {
      H_CI = c(H_CI, quantile(Bootstrap_Pre[,i],quantiles[2]))
      L_CI = c(L_CI, quantile(Bootstrap_Pre[,i],quantiles[1]))
    }
    Boot_highCI[[nn]] = H_CI
    Boot_lowCI[[nn]] = L_CI
  }
  names(Boot_preds) = paste("Bootstrapping ", rownames(bootstrap_results)[1:Result_length],sep = "")
  names(Boot_highCI) = paste("Bootstrapping HighCI ", rownames(bootstrap_results)[1:Result_length],sep = "")
  names(Boot_lowCI) = paste("Bootstrapping LowCI ", rownames(bootstrap_results)[1:Result_length],sep = "")
  attributes = list("Outcome variables" = outcome, "Estimate type" = estimate, "Bootstrapping times" = iterations, "If standardize data" = standardized, "method" = method, "gamma" = gamma, "k" = k, "ktrend" = ktrend  )
  return(list("est" = Boot_preds, "highCI" = Boot_highCI, "lowCI" = Boot_lowCI, "attributes" = attributes))
}










