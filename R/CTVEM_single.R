#' Continuous Time-Varying Effect model for a single-time bootstrapping estimation.
#'
#' This is the CTVEM function which serves for performing bootstrapping estimation. The output of this function is the single estimation of the bootstrapping process, i.e., randomly select data rows from the original dataset, then do estimation.
#' @param data Specify the data frame that contains the interested data, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
#' @param Time The name of the Time column in the data E.G. Time = "Time" (must be specified).
#' @param outcome This is the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2"). If it is NULL, it will consider each variables as outcome once.
#' @param ID The name of the ID column in the data E.G. ID = "ID"
#' @param estimate The relationship which we are interested, estimate = "marginal" or "partial". Default is the "marginal".
#' @param Tpred A vector which indicates that interested time points, e.g. seq(0,30,1)
#' @param quantiles The quantiles to build bootstrapping CI, the default value is c(low_quantile, high_quantile) = c(.025, 0.975)
#' @param boot Indicate if we perform bootstrapping estimation or not. If boot == True, we perfomr bootstrapping estimation. The default value is False
#' @param plot_show The option to suppress the plot outcomes. The default if FALSE which means the plot outcomes will not appear.
#' @param output_type Indicate which output form will be returned. If output_type == "CI", point estimations and corresponding CIs will be returned. If output_type =="PE", only ponit estimation will be returned. The default value is "CI"
#' @param standardized This specifies whether all of the variables (aside from Time) should be standardized. Options are TRUE, FALSE, and "center". TRUE means within-person standardize each variable (aka get the person-centered z-scores), FALSE means use the raw data, "center" means to only within-person mean-center the variables. Default = TRUE. FALSE is not recommended unless you have done these transformations yourself (OPTIONAL)
#' @param method Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".
#' @param gamma This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param k The number of k selection points used in the model for stage 1 (see ?choose.k in mgcv package for more details) The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL, BUT RECOMMENDED)
#' @param ktrend The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)


CTVEM_single <-
  function(data = NULL,
           Time = "Time",
           outcome = NULL,
           ID = "ID",
           estimate = "marginal",
           Tpred = seq(0, 30, 1),
           #datamanipu = "DT",
           plot_show = FALSE,
           boot = FALSE ,
           output_type = "CI",
           standardized = TRUE,
           method = "bam",
           gamma = 1,
           k = 10,
           ktrend = 3,
           quantiles = c(.025, 0.975),
           predictionsend  = NULL) {

  #LOAD NECESSARY PACKAGES
  #library(mgcv) #USED FOR THE PRIMARY ANALYSES
  #library(plyr)
  #library(zoo)
  #library(reshape2)

  # here - make predictions end just slightly longer than Tpred max
  if(is.null(predictionsend)) predictionsend <- max(Tpred)*1.1

  # Check if the data includes the colnames
  if (is.null(colnames(data))){
    stop("ERROR: Please include colnames for your dataset")
  }


  # Get the names of variables
  colnames_data = colnames(data)
  varnames = colnames_data[-c(which(colnames_data == ID), which(colnames_data == Time))]
  if (sum(is.na(unique(data[, ID]))) > 0) {
    # CHECK TO SEE IF THE THERE ARE ANY MISSING ID Variables #
    stop("ERROR: At least one ID is missing. Missing data is not allowed in the ID column. Replace and re-run")
  }
  numberofpeople <- as.numeric(length(unique(data[, ID]))) #NUMBER OF PEOPLE = NUMBER OF UNIQUE IDs in DATAFRAME

  # # Build randomly choosen data
  # Select = sort(sample(seq(1,nrow(data),1),size = nrow(data),replace = T),decreasing = F)
  # data_select = data[Select,]
  # data_select = data.frame(data_select)
  #
  #data[,"Time"]=data[,Time]


  # Prepare some values
  # Notice that the predictions interval is always 1 when we go to do data manipulation and estimations.
  # However, when we return the estimation values from DTVEM, we only care about the interested time points indicated by users
  # predictionstart = Tpred[1]
  # predictionsend = Tpred[length(Tpred)]
  # predictionsinterval = 1
  numberofknots = k
  list_name = c() # Use to give the name of the ouput list

  # Build matrix to store results
  # When the outcomes is not specified, we will consider each variable as outcome once and return all marginal/partial effects. E.g. if there are three variables X1, X2 and X3.
  # the number of total marginal effects we want to consider is 3X3.
  if(is.null(outcome)){
    Result_length = length(varnames)*length(varnames)
  }else{ # with specified outcomes. E.g. three variables X1, X2 and X3. Outcome is c("X1","X2"). We only consider 3X2 = 6 marginal/partial effects
    Result_length = length(outcome)*length(varnames)
  }

  # We just need the point estimation for doing bootstrapping estimation
  #Single_preds = matrix(nrow = Result_length,ncol = length(seq(predictionstart, predictionsend, by = predictionsinterval))) # Contain each marginal effects estimation in order
  Single_preds = vector("list",length = Result_length) # Contain each marginal effects estimation
  Single_highCI = vector("list",length = Result_length) # Contain Upper CI for each marginal effects estimation
  Singel_lowCI = vector("list", length = Result_length) # Contain Lower CI for each marginal effects estimation


  # Estimate marginal effects
  if(estimate == "marginal"){

    # When the outcomes is not specified, we will consider each variable as outcome once and return all marginal effects
    if(is.null(outcome)){
      varnames_mat <- as.matrix(expand.grid(varnames, varnames))
    }else{
      varnames_mat <- as.matrix(expand.grid(varnames, outcome))
    }
    for (i in 1:nrow(varnames_mat)) {
      #i = 1
      input_list = as.list(varnames_mat[i,]) # Get the input list, for the marginal case. There will always be only 2 variables
      differentialtimevaryingpredictors = varnames_mat[i,1] # Take the first variable as predictor
      outcome_mcr = varnames_mat[i,2] # Take the second variable as predictor

      # now do data manipulation
        datamanipulationout = datamanipulation(
          differentialtimevaryingpredictors = differentialtimevaryingpredictors,
          outcome = outcome_mcr,
          data = data,
          ID = ID,
          Time = Time,
          standardized = standardized,
          predictionsend = predictionsend
        )
        lengthcovariates = datamanipulationout$lengthcovariates
        namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
        laglongreducedummy = datamanipulationout$laglongreducedummy


      # Run the CT estimation
      cat(paste("Perform the ",i,"/",nrow(varnames_mat), " time " , estimate , " CTVEM estimation",".\n",sep=""))

      #print(head(laglongreducedummy))
      estout = CTest(
        differentialtimevaryingpredictors = differentialtimevaryingpredictors,
        outcome = outcome_mcr,
        # predictionstart = predictionstart,
        # predictionsend = predictionsend,
        # predictionsinterval = predictionsinterval,
        namesofnewpredictorvariables = namesofnewpredictorvariables,
        laglongreducedummy = laglongreducedummy,
        method = method,
        gamma = gamma,
        numberofknots = numberofknots,
        ktrend = ktrend,
        lengthcovariates = lengthcovariates,
        plot_show = plot_show
      )
      # Do estimation
      # @ Kejin 23 june - check below

      # get model out
      model = estout$mod

      # make dummy frame for getting predictions out of model
       pdat = data.frame(timediff = Tpred,time = 0)
       add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat))
       pdat2 = cbind(pdat, add_mat)
       colnames(pdat2)[3:ncol(pdat2)] = namesofnewpredictorvariables

       # get predictions from model
       predictions = predict(model, pdat2, type = "terms", se = "TRUE")

      # compute CIs
      Single_preds[[i]] = as.vector(predictions$fit[,1])
      Single_highCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[2], lower.tail = T) * as.vector(predictions$se.fit[, 1])
      Singel_lowCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[1], lower.tail = T) *  as.vector(predictions$se.fit[, 1])
      list_name = c(list_name, namesofnewpredictorvariables)
    }
  }

  # @ kejin - make below work the same as above!

  # Estimate partial effects
  if(estimate == "partial"){
    # We can directly apply CTVEM to estimate all partial effects. So we do not need to use expand.grid to separate each effects.
    # Do data manipulation
    if(is.null(outcome)){
      outcome_pcr = varnames
    }else{
      outcome_pcr = outcome
    }
    input_list = as.list(varnames)
    differentialtimevaryingpredictors = varnames

      datamanipulationout = datamanipulation(
        differentialtimevaryingpredictors = differentialtimevaryingpredictors,
        outcome = outcome_pcr,
        data = data,
        ID = ID,
        Time = Time,
        standardized = standardized,
        predictionsend = predictionsend
      )
      lengthcovariates = datamanipulationout$lengthcovariates
      namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
      laglongreducedummy = datamanipulationout$laglongreducedummy


    # Run the CT estimation
    cat(paste("Perform ", estimate , " CTVEM estimation",".\n",sep=""))
    #print(head(laglongreducedummy))
    estout = CTest(
      differentialtimevaryingpredictors = differentialtimevaryingpredictors,
      outcome = outcome_pcr,
      # predictionstart = predictionstart,
      # predictionsend = predictionsend,
      # predictionsinterval = predictionsinterval,
      namesofnewpredictorvariables = namesofnewpredictorvariables,
      laglongreducedummy = laglongreducedummy,
      method = method,
      gamma = gamma,
      numberofknots = numberofknots,
      ktrend = ktrend,
      lengthcovariates = lengthcovariates,
      plot_show = plot_show
    )
    # Do estimation
    # get model out
    model = estout$mod

    # make dummy frame for getting predictions out of model
    pdat = data.frame(timediff = Tpred,time = 0)
    add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat))
    pdat2 = cbind(pdat, add_mat)
    colnames(pdat2)[3:ncol(pdat2)] = namesofnewpredictorvariables

    # get predictions from model
    predictions = predict(model, pdat2, type = "terms", se = "TRUE")


    for (i in 1:length(Single_preds)) {
      Single_preds[[i]] = as.vector(predictions$fit[,i])
      Single_highCI[[i]] = as.vector(predictions$fit[,i])+qnorm(quantiles[2],lower.tail = T)*as.vector(predictions$se.fit[,i])
      Singel_lowCI[[i]] = as.vector(predictions$fit[,i])+qnorm(quantiles[1],lower.tail = T)*as.vector(predictions$se.fit[,i])
    }
    list_name = c(list_name, namesofnewpredictorvariables)
  }

  names(Single_preds) = paste(estimate," ",list_name,sep = "")
  names(Single_highCI) = paste("HighCI ",estimate," " ,list_name,sep = "")
  names(Singel_lowCI) = paste("LowCI ",estimate," ",list_name,sep = "")


  if(boot == TRUE){
    list_names = names(Single_preds)
    returnmatrix = matrix(unlist(Single_preds), nrow = length(Single_preds),byrow = T) # Since we are doing bootstrapping estimation, we only need to return point estimation from the single CTVEM
    rownames(returnmatrix) = list_names
    return(returnmatrix)
  }else{
    attributes = list("Outcome variables" = outcome, "Estimate type" = estimate, "If standardize data" = standardized, "method" = method, "gamma" = gamma, "k" = k, "ktrend" = ktrend  )
    if(output_type == "CI"){
      return(list("est" = Single_preds, "highCI" = Single_highCI, "lowCI" = Singel_lowCI, "attributes" = attributes))
    }else{
      return(list("est" = Single_preds, "attributes" = attributes))
    }

  }

}





