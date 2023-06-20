#' CT estimation of expct
#'
#' PLEASE USE THE expct FUNCTION RATHER THAN THIS UNLESS YOU WOULD LIKE TO SPECIFY THIS MANUALLY.
#' @param differentialtimevaryingpredictors The variables that will be a varying-coefficient of differential time (AKA the lags you want to know what times they predict the outcome). This must be specified as a vector using c("variables here"). e.g. c("X","Y") (REQUIRED)
#' @param outcome This is each of the outcome variables. Specified as outcome="outcomevariablename" for a single variable or outcome=c("outcomevariablename1","outcomevariablename2") (REQUIRED)
#' @param predictionstart The differential time value to start with, default is NULL, and the lowest time difference in the time series will be used (use lower value if you're first value if you're interested in a smaller interval prediction) e.g. predictionstart = 1. If this is not specified and using a continuous time model, make sure to set blockdata = TRUE so that it will be automatically chosen. (OPTIONAL)
#' @param predictionsend The differential time value to end with. This means how long you want your largest time difference in the study to be (i.e. if you wanted to predict up to allow time predictions up to 24 hours and your time intervals were specified in hours, you would set predictionsend = 24). If this is not specified and using a continuous time model, make sure to set blockdata = TRUE so that it will be automatically chosen. (OPTIONAL)
#' @param namesofnewpredictorvariables This is the name of the predictors.
#' @param laglongreducedummy This is the long data output from the data manipulation.
#' @param method Indicate which method will be used to estimate time-varying effetcs. The default value is "bam". Another option is "gam".
#' @param gamma This can be used to change the wiggliness of the model. This can be useful if the model is too smooth (i.e flat). The lower the number the more wiggly this will be (see ?gam in MGCV for more information). The default is equal to 1. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param numberofknots The number of k selection points used in the model for stage 1 (see ?choose.k in mgcv package for more details) (note that this is for the raw data k2 refers to the k for the re-blocked data), default is 10. The ideal k is the maximum number of data points per person, but this slows down DTVEM and is often not required. (OPTIONAL)
#' @param ktrend The number of k selection points used in the model for the time spline (NOTE THAT THIS CONTROLS FOR TIME TRENDS OF THE POPULATION)  (see ?choose.k in mgcv package for more details). Default is 3. (OPTIONAL)
#' @param debug This will print more useless information as it goes along. Only useful for troubleshooting problems. (OPTIONAL, UNCOMMONLY SPECIFIED)
#' @param lengthcovariates the number of covariates (+2)
#' @param plot_show Controls if show the plot
#' @param weights optional vector of weights to be passed to gam/bam
#' @return The output of this function is: The estimated model
#' @export
#'
#'
CTest = function(differentialtimevaryingpredictors = differentialtimevaryingpredictors,
                 outcome = outcome,
                 # predictionstart = predictionstart,
                 # predictionsend = predictionsend,
                 namesofnewpredictorvariables = namesofnewpredictorvariables,
                 laglongreducedummy = laglongreducedummy,
                 method = "bam",
                 gamma = gamma,
                 numberofknots = numberofknots,
                 ktrend = 3,
                 debug = debug,
                 lengthcovariates = lengthcovariates,
                 plot_show = FALSE,
                 weights = NULL) {  #arguments to be fed to mgcv::bam or mgcv::gam
  #library(mgcv)
  #MODEL SPECIFICATION - MULTIVARIATE OUTCOME
  controlvariables = NULL
  predictionsinterval = 1
  for(j in 1:(length(namesofnewpredictorvariables))){
    differentialterm =
      paste("s(timediff, by = ", namesofnewpredictorvariables[j], ", k=", numberofknots, ")", sep ="")
    if(j ==1){
      differentialtermlist = differentialterm
    }else{
      differentialtermlist = paste(differentialtermlist, differentialterm, sep =
                                     " + ")
    }
  }

  #WRITE THE MODEL STATEMENTS - NEW MULTIVARIATE OUTCOME
  if (is.null(controlvariables)==TRUE){ #IF CONTROL VARIABLES ARE NOT SPECIFIED
    # OR edited - check?
    model <-
      as.formula(paste("outcome~", differentialtermlist, ifelse(
        ktrend < 3, "", paste("+s(time,k=", ktrend, ")")
      )))
  }else{ #IF CONTROL VARIABLES ARE  SPECIFIED
    controlvariablelist = paste(controlvariables, collapse = " + ")
    allpredictorvariables = paste(differentialtermlist, controlvariablelist, sep = " + ")
    model <- as.formula(paste("outcome~", allpredictorvariables, "+s(time,k=", ktrend, ")"))
  }
# browser()
  #print(paste("Varying-coefficient model = ",paste(model)))
  #print(model)
  #NOW RUN THE USER SPECIFIED MODEL - MULTIVARIATE
  if(method == "gam") {
    trytest = try(mod <-
                    gam(model,
                        data = laglongreducedummy,
                        na.action = na.omit,
                        gamma = gamma,
                        weights = weights,
                        drop.intercept = TRUE),
                  silent = TRUE)
    if (try(summary(trytest)[2], silent = TRUE)
        == "try-error") {
      stop("Error:CTest Failed. Try decreasing k.")
    }
  } else{
    trytest = try(mod <-
                    bam(
                      model,
                      data = laglongreducedummy,
                      na.action = na.omit,
                      gamma = gamma,
                      discrete = FALSE,
                      method = "fREML",
                      weights = weights,
                      drop.intercept = TRUE,
                    ),
                  silent = TRUE)
    if (try(summary(trytest)[2], silent = TRUE)
        == "try-error") {
      stop("Error:CTest Failed. Try decreasing k.")
    }
  }


  #RANDOM EFFECT DTVEM TO ESTIMATE PERSON-SPECIFIC EFFECTS
  # laglongreducedummy2=laglongreducedummy
  # laglongreducedummy2$ID=factor(as.factor(laglongreducedummy2$ID),ordered =FALSE)
  #
  # out=gam(outcome~ti(timediff,ID,by=X1lagonX1,k=10,bs="fs"),data=laglongreducedummy2[laglongreducedummy2$ID%in%c(1:20),],gamma=1.2)
  # plot(out,select=1,ylim=c(-.2,.4))
  # plot(out,select=2,ylim=c(-.2,.4))
  #
  # longpreds=data.frame(matrix(NA,ncol=2,nrow=10*20))
  # for(i in 1:20){
  #   pdat=data.frame(timediff=1:10,ID=1,X1lagonX1=1)
  #   predictions=predict(out,pdat,type="terms",se=TRUE)
  #   longpreds[(((i-1)*10)+1):(i*10),1]=predictions$fit[,"ti(timediff,ID):X1lagonX1"]
  #   longpreds[(((i-1)*10)+1):(i*10),2]=i
  # }




  #PLOT EACH OF THE VARYING-COEFFICIENT MODELS

  if(plot_show){
    graphnumber=0
    for(i in 1:length(outcome)){
      for(ii in 1:length(differentialtimevaryingpredictors)){
        graphnumber=graphnumber+1
        plot(mod,select=graphnumber,ylab=paste("Beta Coefficient of ",differentialtimevaryingpredictors[ii],"lag on ",outcome[i],sep=""),xlab="Time Differences",main="expct Estimation stage")
      }
    }
  }

  #SET UP PREDICTIONS MATRIX (CALLED pdat)
  #CHANING largestlagdiff to length((predictionstart/predictionsinterval):(predictionsend/predictionsinterval)*predictionsinterval)
  #REVISED PDAT WITH 10x Sampling
  # if (is.null(predictionstart)) {
  #   mindiff <- min(laglongreducedummy$timediff, na.rm = TRUE)
  #   predictionstart <- mindiff
  # }
  # if (is.null(predictionsend)) {
  #   maxdiff <- max(laglongreducedummy$timediff, na.rm = TRUE)
  #   predictionsend <- maxdiff
  # }
  # if (is.null(predictionsinterval)) {
  #   mindiff <- min(laglongreducedummy$timediff, na.rm = TRUE)
  #   #if(isTRUE(blockdata)){ #IF BLOCK DATA IS YES THEN DIVIDE THE MINIMUM DIFFERENCE BY 10 FOR PRECISION
  #   #  predictionsinterval<-mindiff/10
  #   #}else{
  #   predictionsinterval <- mindiff
  #   #}
  # }


  # Tpred = seq(from = predictionstart, to = predictionsend, by = predictionsinterval)

  # #print(paste("long model is =",model))
  # pdat = matrix(NA, nrow = length(Tpred), ncol = lengthcovariates)
  # #print(paste(lengthcovariates,lengthcovariates)
  # #print(paste("pdat ncol=",ncol(pdat)))
  #
  # pdat[,1]=Tpred
  # pdat[,2]=0
  # pdat[,3:(length(differentialtimevaryingpredictors)*length(outcome)+2)]=1
  # if (length(controlvariables) > 0) {
  #   #IF THE COVARIATES ARE INCLUDED, 0 ZERO OUT THE COVARIATES
  #   pdat[, (lengthcovariates - length(controlvariables) + 1):lengthcovariates] = 0
  # }
  # pdat <- data.frame(pdat)
  # names(pdat)[1] <- "timediff"
  # names(pdat)[2] <- "time"
  #
  # tempnumbcount = 2
  #
  # #NAME THE VARIABLES IN PDAT
  # for (j in 1:length(outcome)){
  #   for(jj in 1:length(differentialtimevaryingpredictors)){ #NAME THE DIFFERENTIAL TIME VARYING EFFECT VARIABLES
  #     tempnumbcount=tempnumbcount+1
  #     names(pdat)[tempnumbcount] <- paste(differentialtimevaryingpredictors[jj],"lagon",outcome[j],sep="")
  #   }
  # }
  # tempnumbcount=lengthcovariates-length(controlvariables)*length(outcome)
  # if(length(controlvariables)>0){
  #   for (jjj in 1:length(outcome)){
  #     for (j in 1:length(controlvariables)){ #NAME THE DIFFERENTIAL TIME VARYING EFFECT VARIABLES
  #       tempnumbcount=tempnumbcount+1
  #       jj=j+1+length(differentialtimevaryingpredictors)
  #       names(pdat)[tempnumbcount] <- paste(controlvariables[j],sep="")
  #     }
  #   }
  # }
  #
  #
  # returnlist=list("mod"=mod,"pdat"=pdat)

  return(list("mod"=mod))
}
