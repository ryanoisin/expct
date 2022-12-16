#' Continuous Time-Varying Effect model for a single-time bootstrapping estimation.
#'
#' This is the expct function which serves for performing bootstrapping estimation. The output of this function is the single estimation of the bootstrapping process, i.e., randomly select data rows from the original dataset, then do estimation.
#' @param datset Specify the data frame that contains the interested data, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
#' @param Time The name of the Time column in the data E.G. Time = "Time" (must be specified).
#' @param outcome This is the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2"). If it is NULL, it will consider each variables as outcome once.
#' @param ID The name of the ID column in the data E.G. ID = "ID"
#' @param estimate The relationship which we are interested, estimate = "marginal" or "partial". Default is the "marginal".
#' @param Tpred A vector which indicates that interested time points, e.g. seq(0,30,1)
#' @param quantiles The quantiles to build bootstrapping CI, the default value is c(low_quantile, high_quantile) = c(.025, 0.975)
#' @param boot Indicate if we perform bootstrapping estimation or not. If boot == True, we perfomr bootstrapping estimation. The default value is False
#' @param plot_show The option to suppress the plot outcomes. The default if FALSE which means the plot outcomes will not appear.
#' @param output_type Indicate which output form will be returned. If output_type == "CI", point estimations and corresponding CIs will be returned. If output_type =="PE", only ponit estimation will be returned. If output_type =="SCI", the Simultaneous CIs will be returned. The default value is "CI"
#' @param standardized This specifies whether all of the variables (aside from Time) should be standardized. Options are TRUE, FALSE, and "center". TRUE means within-person standardize each variable (aka get the person-centered z-scores), FALSE means use the raw data, "center" means to only within-person mean-center the variables. Default = TRUE. FALSE is not recommended unless you have done these transformations yourself (OPTIONAL)
#' @param method Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".
#' @param gamma This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param k The number of k selection points used in the model for stage 1 (see ?choose.k in mgcv package for more details) The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL, BUT RECOMMENDED)
#' @param ktrend The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)
#' @param weighting Test option. Calculate weights for stacked data vector?


expct_single <-
  function(datset = NULL,
           Time = "time",
           outcome = NULL,
           ID = "id",
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
           predictionsend  = NULL,
           weighting = FALSE) {

  #LOAD NECESSARY PACKAGES
  #library(mgcv) #USED FOR THE PRIMARY ANALYSES
  #library(plyr)
  #library(zoo)
  #library(reshape2)

  # here - make predictions end just slightly longer than Tpred max
  if(is.null(predictionsend)) predictionsend <- max(Tpred)*1.1

  # Check if the data includes the colnames
  if (is.null(colnames(datset))){
    stop("ERROR: Please include colnames for your dataset")
  }


  # Get the names of variables
  colnames_data = colnames(datset)
  varnames = colnames_data[-c(which(colnames_data == ID), which(colnames_data == Time))]
  if (sum(is.na(unique(datset[, ID]))) > 0) {
    # CHECK TO SEE IF THE THERE ARE ANY MISSING ID Variables #
    stop("ERROR: At least one ID is missing. Missing data is not allowed in the ID column. Replace and re-run")
  }
  numberofpeople <- as.numeric(length(unique(datset[, ID]))) #NUMBER OF PEOPLE = NUMBER OF UNIQUE IDs in DATAFRAME



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
  Single_lowCI = vector("list", length = Result_length) # Contain Lower CI for each marginal effects estimation
  # if(output_type == "LLCI"){
  #   Single_preds_max = vector("list",length = Result_length)
  # }


  # Estimate marginal effects
  if(estimate == "marginal"){

    # When the outcomes is not specified, we will consider each variable as outcome once and return all marginal effects
    if(is.null(outcome)){
      varnames_mat <- as.matrix(expand.grid(varnames, varnames))
    }else{
      varnames_mat <- as.matrix(expand.grid(varnames, outcome))
    }
    # browser()
    for (i in 1:nrow(varnames_mat)) {
      #i = 1
      input_list = as.list(varnames_mat[i,]) # Get the input list, for the marginal case. There will always be only 2 variables
      differentialtimevaryingpredictors = varnames_mat[i,1] # Take the first variable as predictor
      outcome_mcr = varnames_mat[i,2] # Take the second variable as predictor

      # now do data manipulation
        datamanipulationout = datamanipulation(
          differentialtimevaryingpredictors = differentialtimevaryingpredictors,
          outcome = outcome_mcr,
          datset = datset,
          ID = ID,
          Time = Time,
          standardized = standardized,
          predictionsend = predictionsend
        )
        lengthcovariates = datamanipulationout$lengthcovariates
        namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
        laglongreducedummy = datamanipulationout$laglongreducedummy
# browser()

        # oisin added Weights
        if(isTRUE(weighting)){
          datamanipulationout$laglongreducedummy
          timevec <- unique(datamanipulationout$laglongreducedummy$time)
          occurs <- sapply(timevec, function(s) sum(datamanipulationout$laglongreducedummy$time == s))
          weights_tmp <- 1/occurs
          weights <- rep(NA,nrow(datamanipulationout$laglongreducedummy))
          for(qq in 1:length(timevec)){
            weights[
              datamanipulationout$laglongreducedummy$time == timevec[qq]] <- weights_tmp[qq]
          }
        }else{
          weights <- NULL
        }

      # Run the CT estimation
      cat(paste("Perform the ",i,"/",nrow(varnames_mat), " time " , estimate , " expct estimation",".\n",sep=""))

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
        plot_show = plot_show,
        weights = weights
      )


       # get model out
       model = estout$mod

      # make dummy frame for getting predictions out of model
       pdat = data.frame(timediff = Tpred,time = 0)
       add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat))
       pdat2 = cbind(pdat, add_mat)
       colnames(pdat2)[3:ncol(pdat2)] = namesofnewpredictorvariables

       # get predictions from model
       predictions = predict(model, pdat2, type = "terms", se = "TRUE")
       Single_preds[[i]] = as.vector(predictions$fit[,1])
      # compute CIs
      if(output_type == "CI"){
      Single_highCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[2], lower.tail = T) * as.vector(predictions$se.fit[, 1])
      Single_lowCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[1], lower.tail = T) *  as.vector(predictions$se.fit[, 1])
      }

      if(output_type == "SCI"){
        NN = 10000  # The bootstrap times to generate the bias of estimation
        Vb = vcov(model) # Compute the estimated covariance matrix of estimation of parameters
        predictions_ = predict(model, pdat2, se.fit = "TRUE") # get the prediction on the selected grid with fit.se
        se.fit = predictions_$se.fit
        L = mroot(Vb)
        mm = ncol(L)
        mu = rep(0, nrow(Vb))
        BUdiff = mu + L %*% matrix(rnorm(mm*NN), mm, NN) # Resample the value of bias of estimation under normally distributed assumption
        Cg = predict(model, pdat2, type = "lpmatrix")
        simDev = Cg %*% BUdiff
        absDev = abs(sweep(simDev, 1, se.fit, FUN = "/"))
        masd = apply(absDev, 2, max) # Estimate the distribution of the max standardized deviationbetween the true function and the model estimate
        crit = quantile(masd, prob = quantiles[2]) # Compute the critical value for a given confidence level
        Single_highCI[[i]] = as.vector(predictions$fit[,1]) + crit*predictions_$se.fit
        Single_lowCI[[i]] = as.vector(predictions$fit[,1]) - crit*predictions_$se.fit
      }

       # add "large lag SE" @kejin
      # sd_acf = c(sd_acf, sqrt( n^(-1)*(1 + 2*sum(acf[1:ii]^2))) )

      list_name = c(list_name, namesofnewpredictorvariables)
    }

    if(output_type == "LLCI"){ # Compute the se of correlation based on the large_lag method
      Single_preds_all = vector("list",length = length(varnames)*length(varnames)) # compute each single auto correlation effects
      Single_preds_all_names = vector(length = length(varnames)*length(varnames))
      nn = 1
      lag_len = floor(max(Tpred))
      lag_zero = lag_len + 1
      for (ii in c(1:length(varnames))) {
        for (jj in c(1:length(varnames))){
        differentialtimevaryingpredictors = varnames[ii] # Take the first variable as predictor
        outcome_mcr = varnames[jj] # Take the second variable as predictor
        datamanipulationout = datamanipulation(
          differentialtimevaryingpredictors = differentialtimevaryingpredictors,
          outcome = outcome_mcr,
          datset = datset,
          ID = ID,
          Time = Time,
          standardized = standardized,
          predictionsend = predictionsend
        )
        lengthcovariates = datamanipulationout$lengthcovariates
        namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
        laglongreducedummy = datamanipulationout$laglongreducedummy
        estout = CTest(
          differentialtimevaryingpredictors = differentialtimevaryingpredictors,
          outcome = outcome_mcr,
          namesofnewpredictorvariables = namesofnewpredictorvariables,
          laglongreducedummy = laglongreducedummy,
          method = method,
          gamma = gamma,
          numberofknots = numberofknots,
          ktrend = ktrend,
          lengthcovariates = lengthcovariates,
          plot_show = plot_show,
          weights = weights
        )
        model = estout$mod
        Tpred_max = seq(0, max(Tpred),1) # Since user may not be interested in the seq(0,30,1) time points, we estimate correlation effects at seq(0,max(Tpred),1)
        # because the compuation of large_lag se need all previous information
        pdat = data.frame(timediff = Tpred_max,time = 0)
        add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat))
        pdat2 = cbind(pdat, add_mat)
        colnames(pdat2)[3:ncol(pdat2)] = namesofnewpredictorvariables
        # get predictions from model
        predictions_max = predict(model, pdat2, type = "terms", se = "TRUE")
        Single_preds_all[[nn]] = predictions_max$fit[,1]
        Single_preds_all_names[nn] = paste(varnames[ii],"lagon",varnames[jj],sep = "")
        nn = nn + 1
      }
}
      names(Single_preds_all) = Single_preds_all_names
      n = nrow(datset)
      for (iii in c(1:nrow(varnames_mat))) {
        if(length(unique(varnames_mat[iii,])) == 1){ # compute CI of auto correlation effects
          sd_acf = c()
          Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "")
          for(l in Tpred){
            sd_acf = c(sd_acf, sqrt( n^(-1)*(1 + 2*sum(Single_preds_all[[Preds_names]][2:floor(l)]^2))) ) # Compute the large_lag_standard error
          }
          Single_highCI[[iii]]  = Single_preds[[iii]] + 1.96*sd_acf
          Single_lowCI[[iii]] = Single_preds[[iii]] - 1.96*sd_acf
        }else{ # compute CI of cross correlation effects
          c1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "")
          c2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,1],sep = "")
          a1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,1],sep = "")
          a2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,2],sep = "")

          ccf_ij_twoway = c(rev(Single_preds_all[[c2Preds_names]][-1]),Single_preds_all[[c1Preds_names]])
          #ccf_ji_twoway = rev(ccf_ij_twoway)
          acf_i = Single_preds_all[[a1Preds_names]][-1]
          acf_j = Single_preds_all[[a2Preds_names]][-1]
          acf_i_two_way = c(rev(acf_i),1,acf_i)
          acf_j_two_way = c(rev(acf_j),1,acf_j)
          sd_ccf = c()
          for (l in Tpred){
            floor_l = floor(l)
            lag_loc = lag_zero+floor_l
            corr_ij_k_add_v = ccf_ij_twoway[(lag_loc-lag_len) : (lag_loc+lag_len)]
            corr_ij_k_minus_v = ccf_ij_twoway[(lag_loc+lag_len) : (lag_loc-lag_len)]
            #corr_ji_k_add_v = ccf_ji_twoway[(lag_loc-lag_len) : (lag_loc+lag_len)]
            #corr_ji_k_minus_v = ccf_ji_twoway[(lag_loc+lag_len) : (lag_loc-lag_len)]
            #corr_ii_k_add_v = acf_i_two_way[(lag_loc-lag_len) : (lag_loc+lag_len)]
            corr_jj_k_add_v = acf_j_two_way[(lag_loc-lag_len) : (lag_loc+lag_len)]
            #corr_ii_k_add_v[is.na(corr_ii_k_add_v)] = 0
            corr_jj_k_add_v[is.na(corr_jj_k_add_v)] = 0
            corr_ij_k_add_v[is.na(corr_ij_k_add_v)] = 0
            corr_ij_k_minus_v[is.na(corr_ij_k_minus_v)] = 0
            #corr_ji_k_add_v[is.na(corr_ji_k_add_v)] = 0
            #corr_ji_k_minus_v[is.na(corr_ji_k_minus_v)] = 0
            #sd_ccf = c(sd_ccf, sqrt( (n-floor_l)^(-1)*(1 + 2*sum(Single_preds_all[[data_1_name]][1:floor_l]*Single_preds_all[[data_2_name]][1:floor_l])))) # Compute the large_lag_standard error based in Eq 12.1.9 of the book TS analysis and control
            sd_ccf = c(sd_ccf, sqrt((n-floor_l)^(-1)*(
              sum(acf_i_two_way*acf_j_two_way)
              + sum(corr_ij_k_add_v*corr_ij_k_minus_v)
              +  ccf_ij_twoway[lag_loc]^2*(sum(0.5*acf_i_two_way^2)+sum(0.5*acf_j_two_way^2) + sum(ccf_ij_twoway^2))
              -  2*ccf_ij_twoway[lag_loc]*(sum(acf_i_two_way*corr_ij_k_add_v) + sum(rev(ccf_ij_twoway)*corr_jj_k_add_v))
            )))


          }

          Single_highCI[[iii]]= (Single_preds[[iii]] + 1.96*sd_ccf) # Compute the highCI
          Single_lowCI[[iii]]= (Single_preds[[iii]] - 1.96*sd_ccf)
        }

      }
    }

  }

  # Estimate partial effects
  if(estimate == "partial"){
    # We can directly apply expct to estimate all partial effects. So we do not need to use expand.grid to separate each effects.
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
        datset = datset,
        ID = ID,
        Time = Time,
        standardized = standardized,
        predictionsend = predictionsend
      )
      lengthcovariates = datamanipulationout$lengthcovariates
      namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
      laglongreducedummy = datamanipulationout$laglongreducedummy
    # Run the CT estimation
    cat(paste("Perform ", estimate , " expct estimation",".\n",sep=""))
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
    for (i in 1:length(Single_preds)){
      Single_preds[[i]] = as.vector(predictions$fit[,i])
    }
    if(output_type == "CI"){
    for (i in 1:length(Single_preds)){
      Single_highCI[[i]] = as.vector(predictions$fit[,i])+qnorm(quantiles[2],lower.tail = T)*as.vector(predictions$se.fit[,i])
      Single_lowCI[[i]] = as.vector(predictions$fit[,i])+qnorm(quantiles[1],lower.tail = T)*as.vector(predictions$se.fit[,i])
    }
    }

    if(output_type == "SCI"){
      for (i in 1:length(Single_preds)){
        NN = 10000  # The bootstrap times to generate the bias of estimation
        Vb = vcov(model) # Compute the estimated covariance matrix of estimation of parameters
        predictions_ = predict(model, pdat2, se.fit = "TRUE") # get the prediction on the selected grid with fit.se
        se.fit = predictions_$se.fit
        L = mroot(Vb)
        mm = ncol(L)
        mu = rep(0, nrow(Vb))
        BUdiff = mu + L %*% matrix(rnorm(mm*NN), mm, NN) # Resample the value of bias of estimation under normally distributed assumption
        Cg = predict(model, pdat2, type = "lpmatrix")
        simDev = Cg %*% BUdiff
        absDev = abs(sweep(simDev, 1, se.fit, FUN = "/"))
        masd = apply(absDev, 2, max) # Estimate the distribution of the max standardized deviationbetween the true function and the model estimate
        crit = quantile(masd, prob = quantiles[2]) # Compute the critical value for a given confidence level
        Single_highCI[[i]] = as.vector(predictions$fit[,i]) + crit*predictions_$se.fit
        Single_lowCI[[i]] = as.vector(predictions$fit[,i]) - crit*predictions_$se.fit
      }
    }

    list_name = c(list_name, namesofnewpredictorvariables)

    if(output_type == "LLCI"){ # Compute the se of correlation based on the large_lag method partial version
      # We need estimate all partial effects at one time since the large_lag se need these information
      # Get the varmat as we did in the marginal estimation
      if(is.null(outcome)){
        varnames_mat <- as.matrix(expand.grid(varnames, varnames))
      }else{
        varnames_mat <- as.matrix(expand.grid(varnames, outcome))
      }
      pp = length(varnames)
      Single_preds_all = vector("list",length = pp*pp) # record each single auto and cross correlation effects

      #names(Single_preds_all) = varnames
      differentialtimevaryingpredictors = varnames
      datamanipulationout = datamanipulation(
        differentialtimevaryingpredictors = differentialtimevaryingpredictors,
        outcome = varnames,
        datset = datset,
        ID = ID,
        Time = Time,
        standardized = standardized,
        predictionsend = predictionsend
      )
      lengthcovariates = datamanipulationout$lengthcovariates
      namesofnewpredictorvariables = datamanipulationout$namesofnewpredictorvariables
      laglongreducedummy = datamanipulationout$laglongreducedummy
      cat(paste("Perform ", estimate , " expct estimation for large lag CI",".\n",sep=""))
      #print(head(laglongreducedummy))
      estout = CTest(
        differentialtimevaryingpredictors = differentialtimevaryingpredictors,
        outcome = outcome_pcr,
        namesofnewpredictorvariables = namesofnewpredictorvariables,
        laglongreducedummy = laglongreducedummy,
        method = method,
        gamma = gamma,
        numberofknots = numberofknots,
        ktrend = ktrend,
        lengthcovariates = lengthcovariates,
        plot_show = plot_show
      )
      model = estout$mod
      Tpred_max = seq(0, max(Tpred),1)
      pdat = data.frame(timediff = Tpred_max,time = 0)
      add_mat = matrix(1,  ncol = length(namesofnewpredictorvariables),nrow = nrow(pdat))
      pdat2 = cbind(pdat, add_mat)
      colnames(pdat2)[3:ncol(pdat2)] = namesofnewpredictorvariables

      # get predictions from model
      predictions = predict(model, pdat2, type = "terms", se = "TRUE")

      # Store all partial auto and cross correlation
      nn = 1
      Single_preds_all_names = vector(length = pp*pp)
      lag_len = floor(max(Tpred))
      lag_zero = lag_len + 1
      for(i in c(1:pp)){
        for(j in c(1:pp)){
        Single_preds_all[[nn]] = predictions$fit[,paste("s(timediff):",varnames[i],"lagon",varnames[j],sep = "")]
        Single_preds_all_names[nn] = paste(varnames[i],"lagon",varnames[j],sep = "")
        nn = nn + 1
        }
      }
      names(Single_preds_all) = Single_preds_all_names
      n = nrow(datset)
      for (iii in c(1:nrow(varnames_mat))) {
        if(length(unique(varnames_mat[iii,])) == 1){ # compute CI of auto correlation effects
          sd_acf = c()
          Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "")
          for(l in Tpred){
            sd_acf = c(sd_acf, sqrt( n^(-1)*(1 + 2*sum(Single_preds_all[[Preds_names]][2:floor(l)]^2))) ) # Compute the large_lag_standard error
          }
          Single_highCI[[iii]]  = Single_preds[[iii]] + 1.96*sd_acf
          Single_lowCI[[iii]] = Single_preds[[iii]] - 1.96*sd_acf
        }else{ # compute CI of cross correlation effects
          c1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "")
          c2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,1],sep = "")
          a1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,1],sep = "")
          a2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,2],sep = "")

          ccf_ij_twoway = c(rev(Single_preds_all[[c2Preds_names]][-1]),Single_preds_all[[c1Preds_names]])
          #ccf_ji_twoway = rev(ccf_ij_twoway)
          acf_i = Single_preds_all[[a1Preds_names]][-1]
          acf_j = Single_preds_all[[a2Preds_names]][-1]
          acf_i_two_way = c(rev(acf_i),1,acf_i)
          acf_j_two_way = c(rev(acf_j),1,acf_j)
          sd_ccf = c()
          for (l in Tpred){
            floor_l = floor(l)
            lag_loc = lag_zero+floor_l
            corr_ij_k_add_v = ccf_ij_twoway[(lag_loc-lag_len) : (lag_loc+lag_len)]
            corr_ij_k_minus_v = ccf_ij_twoway[(lag_loc+lag_len) : (lag_loc-lag_len)]
            #corr_ji_k_add_v = ccf_ji_twoway[(lag_loc-lag_len) : (lag_loc+lag_len)]
            #corr_ji_k_minus_v = ccf_ji_twoway[(lag_loc+lag_len) : (lag_loc-lag_len)]
            #corr_ii_k_add_v = acf_i_two_way[(lag_loc-lag_len) : (lag_loc+lag_len)]
            corr_jj_k_add_v = acf_j_two_way[(lag_loc-lag_len) : (lag_loc+lag_len)]
            #corr_ii_k_add_v[is.na(corr_ii_k_add_v)] = 0
            corr_jj_k_add_v[is.na(corr_jj_k_add_v)] = 0
            corr_ij_k_add_v[is.na(corr_ij_k_add_v)] = 0
            corr_ij_k_minus_v[is.na(corr_ij_k_minus_v)] = 0
            #corr_ji_k_add_v[is.na(corr_ji_k_add_v)] = 0
            #corr_ji_k_minus_v[is.na(corr_ji_k_minus_v)] = 0
            #sd_ccf = c(sd_ccf, sqrt( (n-floor_l)^(-1)*(1 + 2*sum(Single_preds_all[[data_1_name]][1:floor_l]*Single_preds_all[[data_2_name]][1:floor_l])))) # Compute the large_lag_standard error based in Eq 12.1.9 of the book TS analysis and control
            sd_ccf = c(sd_ccf, sqrt((n-floor_l)^(-1)*(
              sum(acf_i_two_way*acf_j_two_way)
              + sum(corr_ij_k_add_v*corr_ij_k_minus_v)
              +  ccf_ij_twoway[lag_loc]^2*(sum(0.5*acf_i_two_way^2)+sum(0.5*acf_j_two_way^2) + sum(ccf_ij_twoway^2))
              -  2*ccf_ij_twoway[lag_loc]*(sum(acf_i_two_way*corr_ij_k_add_v) + sum(rev(ccf_ij_twoway)*corr_jj_k_add_v))
            )))


          }

          Single_highCI[[iii]]= (Single_preds[[iii]] + 1.96*sd_ccf) # Compute the highCI
          Single_lowCI[[iii]]= (Single_preds[[iii]] - 1.96*sd_ccf)
        }
      }
    }




  }





  names(Single_preds) = paste(estimate," ",list_name,sep = "")
  names(Single_highCI) = paste("HighCI ",estimate," " ,list_name,sep = "")
  names(Single_lowCI) = paste("LowCI ",estimate," ",list_name,sep = "")


  if(boot == TRUE | boot == "MBB"){
    list_names = names(Single_preds)
    returnmatrix = matrix(unlist(Single_preds), nrow = length(Single_preds),byrow = T) # Since we are doing bootstrapping estimation, we only need to return point estimation from the single expct
    rownames(returnmatrix) = list_names
    return(returnmatrix)
  }else{
    attributes = list("Outcome variables" = outcome, "Estimate type" = estimate, "If standardize data" = standardized, "method" = method, "gamma" = gamma, "k" = k, "ktrend" = ktrend  )
    if(output_type == "CI" | output_type == "SCI" | output_type == "LLCI"){
      return(list("est" = Single_preds, "highCI" = Single_highCI, "lowCI" = Single_lowCI,"laglongreducedummy" = laglongreducedummy, "attributes" = attributes))
    }else{
      return(list("est" = Single_preds, "laglongreducedummy" = laglongreducedummy, "attributes" = attributes))
    }

  }

}



#Compute_LLE = function(expct_pre = NULL, varnames_mat = NULL, Single_preds = NULL, Single_highCI    )




