#' Continuous Time-Varying Effect model for a single-time bootstrapping estimation.
#'
#' This is the expct function which serves for performing bootstrapping estimation. The output of this function is the single estimation of the bootstrapping process, i.e., randomly select data rows from the original dataset, then do estimation.
#' @param dataset Specify the data frame that contains the interested data, Time (measuing time) and ID column. MUST INCLUDE COLNAMES.
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
#' @param llc_method determine which large-lag-correlation method used to compute the variance of the auto-correlation and cross-correlation

expct_single <-
  function(dataset = NULL,
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
           weighting = FALSE,
           llc_method = "ks") {

  #LOAD NECESSARY PACKAGES
  #library(mgcv) #USED FOR THE PRIMARY ANALYSES
  #library(plyr)
  #library(zoo)
  #library(reshape2)

  if(output_type == "LLCI"){
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

    # check if integers
    intcheck <- all(sapply(Tpred, is.wholenumber))

    # check if equally spaced
    diffcheck <- all(diff(Tpred) == 1)

    # if either condition fails then
    if(!(intcheck & diffcheck)){

    Tpred <- seq(0, round(max(Tpred),0), 1)

    warning("LLCI computation only implemented for equally spaced integer values of Tpred \\
            Returning results for equally spaced Tpred from 0 to maximum user-specified")
    }
  } #end

  # here - make predictions end just slightly longer than Tpred max
  if(is.null(predictionsend)) predictionsend <- max(Tpred)*1.1

  # Check if the data includes the colnames
  if (is.null(colnames(dataset))){
    stop("ERROR: Please include colnames for your dataset")
  }


  # Get the names of variables
  colnames_data = colnames(dataset)
  varnames = colnames_data[-c(which(colnames_data == ID), which(colnames_data == Time))]
  if (sum(is.na(unique(dataset[, ID]))) > 0) {
    # CHECK TO SEE IF THE THERE ARE ANY MISSING ID Variables #
    stop("ERROR: At least one ID is missing. Missing data is not allowed in the ID column. Replace and re-run")
  }
  numberofpeople <- as.numeric(length(unique(dataset[, ID]))) #NUMBER OF PEOPLE = NUMBER OF UNIQUE IDs in DATAFRAME



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
    Result_length = length(outcome)*length(outcome)
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
      varnames_mat <- as.matrix(expand.grid(outcome, outcome))
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
          dataset = dataset,
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
       predictions = predict(model, pdat2, type = "terms", se.fit = "TRUE")
       Single_preds[[i]] = as.vector(predictions$fit[,1])

       # Transform the regression prediction to the correlation
       if(differentialtimevaryingpredictors != outcome_mcr){
         sd_pre = sd(dataset[,differentialtimevaryingpredictors])
         sd_out = sd(dataset[,outcome_mcr])
         Single_preds[[i]] = Single_preds[[i]]* (sd_pre/sd_out)
       }
      # compute CIs
      if(output_type == "CI"){
      Single_highCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[2], lower.tail = T) * as.vector(predictions$se.fit[, 1])
      Single_lowCI[[i]] = as.vector(predictions$fit[, 1]) + qnorm(quantiles[1], lower.tail = T) *  as.vector(predictions$se.fit[, 1])
      if(differentialtimevaryingpredictors != outcome_mcr){
        sd_pre = sd(dataset[,differentialtimevaryingpredictors])
        sd_out = sd(dataset[,outcome_mcr])
        Single_highCI[[i]] = Single_highCI[[i]]* (sd_pre/sd_out)
        Single_lowCI[[i]] = Single_lowCI[[i]]* (sd_pre/sd_out)
      }

      }

      if(output_type == "SCI"){
        NN = 10000  # The bootstrap times to generate the bias of estimation
        Vb = vcov(model, unconditional = TRUE) # Compute the estimated covariance matrix of estimation of parameters
        predictions_ = predict(model, pdat2,type = "terms", se.fit = "TRUE") # get the prediction on the selected grid with fit.se
        se.fit = predictions_$se.fit
        #L = mroot(Vb)
        #mm = ncol(L)
        #mu = rep(0, nrow(Vb))
        #BUdiff = mu + L %*% matrix(rnorm(mm*NN), mm, NN) # Resample the value of bias of estimation under normally distributed assumption

        BUdiff = rmvn(NN,mu = rep(0,nrow(Vb)), V = Vb)

        Cg = predict(model, pdat2, type = "lpmatrix")
        simDev = Cg %*% t(BUdiff)
        absDev = abs(sweep(simDev, 1, se.fit[,1], FUN = "/"))
        masd = apply(absDev, 2, max) # Estimate the distribution of the max standardized deviationbetween the true function and the model estimate
        crit_up = quantile(masd, prob = quantiles[2]) # Compute the critical value for a given confidence level


        Single_highCI[[i]] = as.vector(predictions$fit[,1]) + crit_up*predictions_$se.fit
        Single_lowCI[[i]] = as.vector(predictions$fit[,1]) - crit_up*predictions_$se.fit
        if(differentialtimevaryingpredictors != outcome_mcr){
          sd_pre = sd(dataset[,differentialtimevaryingpredictors])
          sd_out = sd(dataset[,outcome_mcr])
          Single_highCI[[i]] = Single_highCI[[i]]* (sd_pre/sd_out)
          Single_lowCI[[i]] = Single_lowCI[[i]]* (sd_pre/sd_out)
        }
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
        outcome_mcr = varnames[jj] # Take the second variable as response
        datamanipulationout = datamanipulation(
          differentialtimevaryingpredictors = differentialtimevaryingpredictors,
          outcome = outcome_mcr,
          dataset = dataset,
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
      n = nrow(dataset)
      # make est matrix which contains need rx ry rxy and lag order
      for (iii in c(1:nrow(varnames_mat))) {
        if(length(unique(varnames_mat[iii,])) == 1){
          Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "") # get the prediction name
          acf_x = Single_preds_all[[Preds_names]]
          acf_x2 = c(rev(acf_x[-1]),1,acf_x[-1])
          lagx2 = seq(-max(Tpred),max(Tpred),1)
          ccfmat2 = cbind(acf_x2, lagx2)
          ests2 = as.data.frame(cbind(acf_x2, acf_x2, ccfmat2))
          colnames(ests2) = c("rx", "ry", "rxy", "lag")

          # Create the row with NA lag and rx = ry = rxy = 0
          ests2 = rbind(ests2, 0)
          ests2[nrow(ests2),ncol(ests2)] = NA

          # now calculate this over the vector Tpred

          sdvec <- sapply(Tpred, function(k){
            sigma_initial = 1
            lag.max = max(Tpred)
            ivec <- seq(-lag.max, lag.max, 1)
            varests = 0
            itr = 0
            if(llc_method == "kc"){
              while (varests<=0 & itr < 50) {
                varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_kc(ests = ests2, i = s,k, lag_max = lag.max, sigma_sq = sigma_initial)))
                sigma_initial = sigma_initial/2
                itr = itr + 1
              }

            }
            if(llc_method == "ks"){
              while (varests<=0 & itr < 50) {
                varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_ks(ests = ests2, i = s,k, lag_max = lag.max, sigma_sq = sigma_initial)))
                sigma_initial = sigma_initial/2
                itr = itr + 1
              }
            }

            if(llc_method == "c"){
              varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_c(ests = ests2, i = s,k)))
            }

            if(llc_method == "s"){
              varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_s(ests = ests2, i = s,k)))
            }

            sqrt(varests)
          }
          )

          Single_highCI[[iii]] = Single_preds[[iii]] + 1.96*sdvec
          Single_lowCI[[iii]] = Single_preds[[iii]] - 1.96*sdvec

        }else{
          c1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,2],sep = "") # Here location 1 represents the predictor
          c2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,1],sep = "")
          a1Preds_names = paste(varnames_mat[iii,1],"lagon",varnames_mat[iii,1],sep = "")
          a2Preds_names = paste(varnames_mat[iii,2],"lagon",varnames_mat[iii,2],sep = "")

          sigma_pre = sd(dataset[,varnames_mat[iii,1]])
          sigma_out = sd(dataset[,varnames_mat[iii,2]])
          cor_pre_out = cor(dataset[,varnames_mat[iii,1]],dataset[,varnames_mat[iii,2]])

          cr_est = Single_preds_all[[c1Preds_names]]
          ccf_est = cr_est[-1]*(sigma_pre/sigma_out)

          rcr_est = Single_preds_all[[c2Preds_names]]
          rccf_est = rev(rcr_est[-1]*(sigma_out / sigma_pre))
          crossmat = cbind(c(rccf_est, cor_pre_out, ccf_est),seq(-max(Tpred),max(Tpred),1))  # Can Tpred be not an integer;

          acf_pre = Single_preds_all[[a1Preds_names]]
          acf_out = Single_preds_all[[a2Preds_names]]

          ests <- as.data.frame(cbind(c(rev(acf_pre[-1]),1,acf_pre[-1]), c(rev(acf_out[-1]),1,acf_out[-1]), crossmat))
          colnames(ests) <- c("rx", "ry", "rxy", "lag")

          # Create the row with NA lag and rx = ry = rxy = 0
          ests = rbind(ests, 0)
          ests[nrow(ests),ncol(ests)] = NA


          # now calculate this over a vector of k values

          sdvec <- sapply(Tpred, function(k){
            sigma_initial = 1
            lag.max = max(Tpred)
            ivec <- seq(-lag.max, lag.max, 1)
            varests = 0
            itr = 0
            if(llc_method == "kc"){
            while (varests<=0 & itr < 50) {
              varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_kc(ests = ests, i = s,k, lag_max = lag.max, sigma_sq = sigma_initial)))
              sigma_initial = sigma_initial/2
              itr = itr + 1
            }

            }
            if(llc_method == "ks"){
              while (varests<=0 & itr < 50) {
                varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_ks(ests = ests, i = s,k, lag_max = lag.max, sigma_sq = sigma_initial)))
                sigma_initial = sigma_initial/2
                itr = itr + 1
              }
            }

            if(llc_method == "c"){
              varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_c(ests = ests, i = s,k)))
            }

            if(llc_method == "s"){
              varests = (1/(n-k))*sum(sapply(ivec, function(s) var_piece_s(ests = ests, i = s,k)))
            }

            sqrt(varests)
            }
            )



          Single_highCI[[iii]] = Single_preds[[iii]] + 1.96*sdvec
          Single_lowCI[[iii]] = Single_preds[[iii]] - 1.96*sdvec
        }

      }
    }

  }

  # Estimate partial effects
  if(estimate == "partial"){
    warning("partial correlation effects estimation has not been considered so far")
    break
 
}


# Rename the output name from the format Y1lagonY2 to Y1toY2, here Y1 is the predictor and Y2 is the response variable

  list_name_new = gsub("lagon", "to", list_name)

  names(Single_preds) = list_name_new
  names(Single_highCI) = list_name_new
  names(Single_lowCI) = list_name_new


  if(boot == TRUE | boot == "MBB"){
    list_names = names(Single_preds)
    returnmatrix = matrix(unlist(Single_preds), nrow = length(Single_preds),byrow = T) # Since we are doing bootstrapping estimation, we only need to return point estimation from the single expct
    rownames(returnmatrix) = list_names
    return(returnmatrix)
  }else{
    attributes = list("Outcome variables" = outcome, "Estimate type" = estimate, "If standardize data" = standardized, "method" = method, "gamma" = gamma, "k" = k, "ktrend" = ktrend, "Tpred" = Tpred  )
    if(output_type == "CI" | output_type == "SCI" | output_type == "LLCI"){
      return(list("est" = Single_preds, "highCI" = Single_highCI, "lowCI" = Single_lowCI,"laglongreducedummy" = laglongreducedummy, "attributes" = attributes))
    }else{
      return(list("est" = Single_preds, "laglongreducedummy" = laglongreducedummy, "attributes" = attributes))
    }

  }

}



#Compute_LLE = function(expct_pre = NULL, varnames_mat = NULL, Single_preds = NULL, Single_highCI    )




var_piece_c <- function(ests,i,k){

  # figure out what row of the estimates matrix to use for what index
  ir <- which(ests$lag == i) ; if(length(ir)==0){  ir <- which(is.na(ests$lag)) }
  imin <- which(ests$lag == -i) ;  if(length(imin)==0){  imin <- which(is.na(ests$lag)) }
  kr <- which(ests$lag == k); if(length(kr)==0){  kr <- which(is.na(ests$lag)) }
  i_min_k <- which(ests$lag == -i + k); if(length(i_min_k)==0){  i_min_k <- which(is.na(ests$lag)) }
  i_plus_k <- which(ests$lag == i + k); if(length(i_plus_k)==0){  i_plus_k <- which(is.na(ests$lag)) }

  # evaluate expression

  {ests$rx[ir]*ests$ry[ir] + ests$rxy[i_min_k]*ests$rxy[i_plus_k]  -
      2*ests$rxy[kr]*(ests$rx[ir]*ests$rxy[i_plus_k] + ests$rxy[imin]*ests$ry[i_plus_k]) +
      (ests$rxy[kr]^2)*(ests$rxy[ir]^2 + 0.5*((ests$rx[ir])^2) + 0.5*((ests$ry[ir])^2))
  }
}

var_piece_s <- function(ests,i,k){

  # figure out what row of the estimates matrix to use for what index
  ir <- which(ests$lag == i) ; if(length(ir)==0){  ir <- which(is.na(ests$lag)) }
  imin <- which(ests$lag == -i) ;  if(length(imin)==0){  imin <- which(is.na(ests$lag)) }
  kr <- which(ests$lag == k); if(length(kr)==0){  kr <- which(is.na(ests$lag)) }
  i_min_k <- which(ests$lag == -i + k); if(length(i_min_k)==0){  i_min_k <- which(is.na(ests$lag)) }
  i_plus_k <- which(ests$lag == i + k); if(length(i_plus_k)==0){  i_plus_k <- which(is.na(ests$lag)) }

  # evaluate expression

  {ests$rx[ir]*ests$ry[ir]
  }
}


var_piece_ks <- function(ests,i,k,lag_max,sigma_sq){ # kernel estimator with simple expression

  # figure out what row of the estimates matrix to use for what index
  ir <- which(ests$lag == i) ; if(length(ir)==0){  ir <- which(is.na(ests$lag)) }
  imin <- which(ests$lag == -i) ;  if(length(imin)==0){  imin <- which(is.na(ests$lag)) }
  kr <- which(ests$lag == k); if(length(kr)==0){  kr <- which(is.na(ests$lag)) }
  i_min_k <- which(ests$lag == -i + k); if(length(i_min_k)==0){  i_min_k <- which(is.na(ests$lag)) }
  i_plus_k <- which(ests$lag == i + k); if(length(i_plus_k)==0){  i_plus_k <- which(is.na(ests$lag)) }

  # evaluate expression

  {exp(-(i/lag_max)^2/sigma_sq)*ests$rx[ir]*ests$ry[ir]}
}

var_piece_kc <- function(ests,i,k,lag_max,sigma_sq){ # kernel estimator with complete expression

  # figure out what row of the estimates matrix to use for what index
  ir <- which(ests$lag == i) ; if(length(ir)==0){  ir <- which(is.na(ests$lag)) }
  imin <- which(ests$lag == -i) ;  if(length(imin)==0){  imin <- which(is.na(ests$lag)) }
  kr <- which(ests$lag == k); if(length(kr)==0){  kr <- which(is.na(ests$lag)) }
  i_min_k <- which(ests$lag == -i + k); if(length(i_min_k)==0){  i_min_k <- which(is.na(ests$lag)) }
  i_plus_k <- which(ests$lag == i + k); if(length(i_plus_k)==0){  i_plus_k <- which(is.na(ests$lag)) }

  # evaluate expression

  { exp(-(i/lag_max)^2/sigma_sq)*(ests$rx[ir]*ests$ry[ir] + ests$rxy[i_min_k]*ests$rxy[i_plus_k]  - 2*ests$rxy[kr]*(ests$rx[ir]*ests$rxy[i_plus_k] + ests$rxy[imin]*ests$ry[i_plus_k]) +(ests$rxy[kr]^2)*(ests$rxy[ir]^2 + 0.5*((ests$rx[ir])^2) + 0.5*((ests$ry[ir])^2)))}
}


