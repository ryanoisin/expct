# to dos:
# - suppress detrending? k3 = FALSE? k3 = 0? check CTest (line 50) %%Solution: by adding an argument time_trend, if time_trend = FALSE, k3 = 0 %%FINISHED
# - make evertyhing into one function call %%Solution: add a boot argument in the CT_LAG function.
# - make version control on github  %%FINISHED
# - output looking the same %% FINISHED
# - efficiency (nick suggestions) %% Add an argument bam to indicate if we use the bam function to run these code.%% FINISHED
# - test bootstrap feasibility



#' Continuous Time-Varying Effect model (CTVEM)
#'
#' This is the PRIMARY CTVEM function which works for providing some lag information about continuous data.
#' @param data Specify the data frame that contains the interested variables, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
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
#' @param datamanipu (can delete later)Determins which data manipulation method will be used. "DT" means the old one; Otherwise, apply the new one (for the moment, works for only 1 person case)
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
#' @examples
#'

CTVEM <- function(data = NULL,
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
                  weighting = FALSE

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
    Result = CTVEM_single(
      data = data,
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
      weighting = weighting
    )
  }else if(boot == TRUE | boot  == "MBB" ){
    Result = CTVEM_boot(
      data = data,
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
      ctype = ctype
    )
  }

  return(Result)
  #
  # #data[,"Time"]=data[,Time]
  #
  # # Build lists to store results
  #   # When the outcomes is not specified, we will consider each variable as outcome once and return all marginal/partial effects. E.g. if there are three variables X1, X2 and X3.
  #   # the number of total marginal effects we want to consider is 3X3.
  #   if(is.null(outcome)){
  #     Result_length = length(varnames)*length(varnames)
  #   }else{ # with specified outcomes. E.g. three variables X1, X2 and X3. Outcome is c("X1","X2"). We only consider 3X2 = 6 marginal/partial effects
  #     Result_length = length(outcome)*length(varnames)
  #   }
  #
  #  Single_preds = vector("list",length = Result_length) # Contain each marginal effects estimation
  #  Single_highCI = vector("list",length = Result_length) # Contain Upper CI for each marginal effects estimation
  #  Singel_lowCI = vector("list", length = Result_length) # Contain Lower CI for each marginal effects estimation
  #
  # # Prepare some values
  #  predictionstart = Tpred[1]
  #  predictionsend = Tpred[2]
  #  predictionsinterval = Tpred[3]
  #  numberofknots = k
  #  list_name = c() # Use to give the name of the ouput list
  #
  # # Estimate marginal effects
  # if(estimate == "marginal"){
  #
  #   # When the outcomes is not specified, we will consider each variable as outcome once and return all marginal effects
  #   if(is.null(outcome)){
  #     varnames_mat <- as.matrix(expand.grid(varnames, varnames))
  #   }else{
  #     varnames_mat <- as.matrix(expand.grid(varnames, outcome))
  #   }
  #     for (i in 1:nrow(varnames_mat)) {
  #       #i = 1
  #       input_list = as.list(varnames_mat[i,]) # Get the input list, for the marginal case. There will always be only 2 variables
  #       differentialtimevaryingpredictors = varnames_mat[i,1] # Take the first variable as predictor
  #       outcome_mcr = varnames_mat[i,2] # Take the second variable as predictor
  #       # Do data manipulation
  #       datamanipulationout = datamanipulation(input_list=input_list,differentialtimevaryingpredictors=differentialtimevaryingpredictors,outcome=outcome_mcr,controlvariables=controlvariables,data=data,ID=ID,Time=Time,controllag=controllag,standardized=standardized,predictionstart=predictionstart,predictionsend=predictionsend,predictionsinterval=predictionsinterval)
  #       lengthcovariates = datamanipulationout$lengthcovariates
  #       namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
  #       laglongreducedummy = datamanipulationout$laglongreducedummy
  #       # Run the CT estimation
  #       cat(paste("Perform ", estimate , " CTVEM estimation",".\n",sep=""))
  #       estout = CTest(
  #         differentialtimevaryingpredictors = differentialtimevaryingpredictors,
  #         outcome = outcome_mcr,
  #         predictionstart = predictionstart,
  #         predictionsend = predictionsend,
  #         predictionsinterval = predictionsinterval,
  #         namesofnewpredictorvariables = namesofnewpredictorvariables,
  #         laglongreducedummy = laglongreducedummy,
  #         gamma = gamma,
  #         numberofknots = numberofknots,
  #         k3 = k3,
  #         controlvariables = controlvariables,
  #         debug = debug,
  #         lengthcovariates = lengthcovariates,
  #         plot_show = plot_show
  #       )
  #       # Do estimation
  #       model = estout$mod
  #       pdat = estout$pdat
  #       pdat2 = data.frame(timediff = seq(predictionstart, predictionsend, by = predictionsinterval),time = 0)
  #       add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat2))
  #       pdat3 = cbind(pdat2, add_mat)
  #       colnames(pdat3)[3:ncol(pdat3)] = namesofnewpredictorvariables
  #       predictions=predict(model,pdat3,type="terms",se="TRUE")
  #       Single_preds[[i]] = as.vector(predictions$fit[,1])
  #       Single_highCI[[i]] = as.vector(predictions$fit[,1])+qnorm(signif_level/2,lower.tail = F)*as.vector(predictions$se.fit[,1])
  #       Singel_lowCI[[i]] = as.vector(predictions$fit[,1])+qnorm(signif_level/2,lower.tail = T)*as.vector(predictions$se.fit[,1])
  #       list_name = c(list_name, namesofnewpredictorvariables)
  #     }
  # }
  #
  #  # Estimate partial effects
  #  if(estimate == "partial"){
  #    # We can directly apply CTVEM to estimate all partial effects. So we do not need to use expand.grid to separate each effects.
  #    # Do data manipulation
  #    if(is.null(outcome)){
  #      outcome_pcr = varnames
  #    }else{
  #      outcome_pcr = outcome
  #    }
  #    input_list = as.list(varnames)
  #    differentialtimevaryingpredictors = varnames
  #    datamanipulationout = datamanipulation(input_list=input_list,differentialtimevaryingpredictors=differentialtimevaryingpredictors,outcome=outcome_pcr,controlvariables=controlvariables,data=data,ID=ID,Time=Time,controllag=controllag,standardized=standardized,predictionstart=predictionstart,predictionsend=predictionsend,predictionsinterval=predictionsinterval)
  #    lengthcovariates = datamanipulationout$lengthcovariates
  #    namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
  #    laglongreducedummy = datamanipulationout$laglongreducedummy
  #    # Run the CT estimation
  #    cat(paste("Perform ", estimate , " CTVEM estimation",".\n",sep=""))
  #    estout = CTest(differentialtimevaryingpredictors=differentialtimevaryingpredictors,outcome=outcome_pcr,predictionstart=predictionstart,predictionsend=predictionsend,predictionsinterval=predictionsinterval,namesofnewpredictorvariables=namesofnewpredictorvariables,laglongreducedummy=laglongreducedummy,gamma=gamma,numberofknots=numberofknots,k3=k3,controlvariables=controlvariables,debug=debug,lengthcovariates=lengthcovariates,plot_show = plot_show)
  #    # Do estimation
  #    model = estout$mod
  #    pdat = estout$pdat
  #    pdat2 = data.frame(timediff = seq(0, predictionsend, by = predictionsinterval),time = 0)
  #    add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat2))
  #    pdat3 = cbind(pdat2, add_mat)
  #    colnames(pdat3)[3:ncol(pdat3)] = namesofnewpredictorvariables
  #    predictions=predict(model,pdat3,type="terms",se="TRUE")
  #    for (i in 1:length(Single_preds)) {
  #    Single_preds[[i]] = as.vector(predictions$fit[,i])
  #    Single_highCI[[i]] = as.vector(predictions$fit[,i])+qnorm(signif_level/2,lower.tail = F)*as.vector(predictions$se.fit[,i])
  #    Singel_lowCI[[i]] = as.vector(predictions$fit[,i])+qnorm(signif_level/2,lower.tail = T)*as.vector(predictions$se.fit[,i])
  #    }
  #    list_name = c(list_name, namesofnewpredictorvariables)
  #  }
  #
  #  names(Single_preds) = paste(estimate,"-",list_name,sep = "")
  #  names(Single_highCI) = paste(estimate,"-HighCI-",list_name,sep = "")
  #  names(Singel_lowCI) = paste(estimate,"-LowCI-",list_name,sep = "")
  # if(SE){
  #   returnlist=list("est" = Single_preds, "highCI" = Single_highCI, "lowCI" = Singel_lowCI)
  # }else{
  #   returnlist=list("est" = Single_preds)
  # }
  # #class(returnlist)<-"DTVEM"
  # return(returnlist)
}




