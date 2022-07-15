#' Data manipulation
#'
#' PLEASE USE THE LAG FUNCTION RATHER THAN THIS UNLESS YOU WOULD LIKE TO SPECIFY THIS MANUALLY. Manipulates the data so that it creates a vector that is both in wide and long formats
#' @param input_list identical to the elipses above, A list of variable names used in the function e.g. "X","Y"
#' @param differentialtimevaryingpredictors A list of predictor variables to be used as predictor variables e.g. c("X","Y")
#' @param outcome The outcome(s) to be predicted e.g. c("X","Y")
#' @param data The data frame to be converted
#' @param ID The name of the ID variable. E.G. ID = "ID" (must be specified).
#' @param Time The name of the Time variable. E.G. Time = "Time" (must be specified).
#' @param standardized This specifies whether all of the variables (aside from Time) should be standardized. Options are TRUE, FALSE, and "center". TRUE means within-person standardize each variable, FALSE means use the raw data, "center" means to only within-person center the variables. Default = TRUE
#' @param predictionstart The differential time value to start with, default is NULL, and the lowest time difference in the time series will be used (use lower value if you're first value if you're interested in a smaller interval prediction)
#' @param predictionsend The differential time value to end with. This should usually be set by the user.
#' @return The output of this will be (1) The data in a very long stacked format (called laglongreducedummy), (2) the data in a wide format with all lags for all variables, time, and time differences (called Timelagsdummy), (3) the name of the differentialtimevarying predictors (namesofnewpredictorvariables), and (4) the length of the long matrix (laglongmatrixlength)
#' @export

datamanipulation_old <- function(input_list = input_list,
                            differentialtimevaryingpredictors = NULL,
                            outcome = NULL,
                            data = NULL,
                            ID = "ID",
                            Time = NULL,
                            standardized = TRUE,
                            predictionstart = NULL,
                            predictionsend = NULL) {
  #library(plyr)
  #library(zoo)
  #library(reshape2)
  predictionsinterval = 1
  controlvariables = NULL
  controllag = NULL
  numberofvars = length(input_list)
  #data<-data[with(data, order(ID, Time)), ]
  newdata <- matrix(NA, nrow = nrow(data), ncol = (numberofvars + 2)) #CREATING THIS FOR THE STANDARDIZED data

  #CREATE THE ID & TIME COLUMN IN A NEW DATASET
  newdata[,1]<-data[,ID] #ID Column
  newdata[,2]<-data[,Time] #Time Column

  #NAME ID & TIME COLUMNS
  names(newdata)[1] <- "ID"
  names(newdata)[2] <- "Time"

  #if(is.null(predictionstart)==FALSE&is.null(predictionsend)==FALSE){
  #  maxobs=length()
  #}

  #COUNT THE NUMBER OF ROWS BY ID
  df_counts <- melt(table(newdata[,1]))  #NUMBER OF ROWS BY ID
  maxobs=max(df_counts[,2]) #Get the max amount of rows for each lag
  #maxobs=NROW(unique(data[,Time])) #try to define the max lags as the maximum possible number of unique rows in the time variable
  maxobsplusone=maxobs+1
  minlags<-min(df_counts[,2])-1
  maxlags<-maxobs-1 #Max number of lags possible
  negmaxlags=maxlags*-1 #Number of max lags negative
  largestlagdiff=max(data[,Time],na.rm=TRUE)-min(data[,Time],na.rm=TRUE) #GET THE MAXIMUM DIFFERENCES OF THE LARGEST LAG AND THE SMALLEST LAG
  #print(largestlagdiff)


  #if(maxlags < numberofknots) {
  #  numberofknots <- maxlags
  #}

  #CREATE A data FRAME FOR NEWDATA
  newdata2<-data.frame(newdata) #NEWDATA AND NEWDATA ARE BASICALLY TEMPORARY MATRICES THAT JUST STORE ID AND TIME FOR THE Timelags data Frame to FOllow
  #CREATE TIME VARIABLE TIME LAGS


  #THIS CREATE TIME LAGS BY EACH ID
  Timelags <- ddply(.data = newdata2, .variables = .(newdata[,1]), function(x){
    data.frame(lag(zoo(x[,2]), k = 0:negmaxlags))}) #2014.07.24: EDITING BECAUSE OF ERRRORS -- REMOVING SORT

  #FOR LONG VERSION OF LAG
  laglongunitlength=nrow(Timelags) #laglongunitlength = The number of rows in Timelags
  laglongmatrixlength=laglongunitlength*maxobs #laglongmatrixlength = NUMBER OF ROWS x NUMBER OF OBSERVATIONS
  if (!is.null(controllag)){
    lagset=length(controllag) #COUNT THE NUMBER OF CONTROLS
  } else{
    lagset=0
  }
  laglongcol=3+numberofvars*2+lagset #COLUMNS: 1-3: (1) ID, (2) Time, (3) Time Differences, 4 - 5...x2 Predictors and Outcome with all Lags, (Last Columns) Control variables
  laglong<-matrix(NA,nrow=laglongmatrixlength,ncol=laglongcol) #CREATE THE MATRIX

  #COUNT TOTAL NUMBER OF PREDICTORS IN MODEL FOR PREDICTIONS
  lengthcovariates<-length(differentialtimevaryingpredictors)*length(outcome)+length(controlvariables)*length(outcome)+2 #+1 IS BECAUSE OF TIMEDIFF

  #if(is.null(predictionsend)){ #IF THE PREDICTIONS END POINT IS NOT SET, SET IT AT THE LARGEST POSSIBLE LAG DIFFERENCE
  #  predictionsend=largestlagdiff
  #}





  #CREATE TIME LAG DIFFERENCES
  for (i in 2:maxobsplusone) #FOR the second column to the max number of lag columns
  {
    im = i+1 #ONE LAG LATER
    ipm = i+maxobs #maxobs because we start with i being 2 because the id variable is 1
    Timelags[,ipm]<-Timelags[,2]-Timelags[,i]#EXPERIMENTING WITH I INSTEAD OF IM
  }

  #CREATE NA COLUMNS FOR CONTROLLAGVARIABLE


  ipmm=ipm+1

  for (i in 1:numberofvars)
  {
    cat(paste("Setting up the data for variable #",i," of ",numberofvars,".\n",sep=""))
    variablename<-unlist(input_list[i])
    variablename2<-as.name(variablename)
    datacolumn<-data[,variablename]
    idcolumn<-data[,ID]
    #print(data[1,testing])
    columnnumber<-as.numeric(match(variablename,names(data)))
    IDcolumnnumber<-as.numeric(match(ID,names(data)))
    #print(columnnumber)
    lastcolumn<-as.numeric(ncol(data))+1
    #STANDARDIZE THE data
    if(standardized==TRUE){
      standardizeddata<-ave(datacolumn, idcolumn, FUN = function(x){scale(x,scale=TRUE)}) #STANDARDIZING THE data BY ID
    } else if(standardized=="center"){
      standardizeddata<-ave(datacolumn, idcolumn, FUN = function(x){scale(x,scale=FALSE)}) #STANDARDIZING THE data BY ID
    } else{
      standardizeddata<-datacolumn
    }
    #print(Time)
    j = i +2
    newdata[,j]<-standardizeddata #ADDING IT TO A NEW data MATRIX
    names(newdata)[j] <- variablename

    # CREATE data FRAME FOR NEWDATA
    newdata2 = data.frame(newdata)
    # LEFT OFF HERE NOT WORKING RIGHT FOR IPMM
    #colnames(Timelags) #CHECK SUCCESS WITH THIS
    ipmm = ipm + 1 + (i - 1) * maxobs #THIS IS TO CALCULATE THE BEGINNING OF THE NEXT MATRIX SERIES
    ipmmm = ipmm + maxobs
    Timelags[,ipmm:ipmmm] <- ddply(.data = newdata2, .variables = .(newdata[,1]), function(x){
      data.frame(lag(zoo(x[,j]), k = 0:negmaxlags))})

    #JULY 6 2016 CODE
    #NEED TO ADD RENAMING HERE
    #names(Timelags)[ipmm:(ipmmm-1)]=(0:maxlags)*predictionsinterval
    #/END JULY 6 2016 CODE

    Timelags[ipmm] <- NULL #DROP THE ID VARIABLE AS ITS REPEATED



    ipmmmm = ipmm - 1
    ipmmmmmm = ipmm + maxlags#+length(((predictionstart/predictionsinterval):((predictionsend)/predictionsinterval)*predictionsinterval))

    #JULY 6 2016 CODE
    #NEED TO ADD RENAMING CONVENTION HERE
    #OLD CODE: names(Timelags)[ipmm:ipmmmmmm]<- paste(variablename2,"lag",c(0:maxlags),sep="")#RENAME ALL VARIABLES SO THAT THE NAMES ARE BETTER DESCRIPTIONS OF THE data

    #if(beforeblock){ #THIS NEW CODE JULY 11, 2016
    #  names(Timelags)[ipmm:ipmmmmmm]<- paste(variablename2,"lag",c(0:maxlags),sep="") #OLD CODE
    #} else{#NEW CODE
    names(Timelags)[ipmm:ipmmmmmm] <-
      paste(variablename2, "lag", (c(0:maxlags) * predictionsinterval), sep =  "")#THIS NEW CODE JULY 11, 2016
    #}#NEW CODE

    #NEW CODE, MAY CAUSE ERROR: names(Timelags)[ipmm:ipmmmmmm]<- paste(variablename2,"lag",(c(0:maxlags)*predictionsinterval),sep="")#RENAME ALL VARIABLES SO THAT THE NAMES ARE BETTER DESCRIPTIONS OF THE data

    names(Timelags)[1]<-"ID"
    names(Timelags)[2]<- "time"

    ####ATTEMPTING WIDE TO LONG FORMAT####
    for (p in 1:maxobs){
      inext=laglongunitlength*(p-1)+1
      inextnext=(p-1)*laglongunitlength+laglongunitlength##LEFT OFF HERE
      ip=p+1
      im=p+maxobs+1
      imm=p+1+(1+i)*maxobs
      if (!is.null(controllag)){
        immm=imm+controllag #CONTROL LAG CHANGE
      }
      stablecol=(1+i)*maxobs+2
      newstablecol=2*(i-1)+4
      newvariablecol=2*(i-1)+5
      laglong[inext:inextnext,1]<-Timelags[,1] #THIS IS FOR THE ID
      laglong[inext:inextnext,2]<-Timelags[,2] #THIS IS FOR THE TIME
      laglong[inext:inextnext,3]<-Timelags[,im] #THIS IS FOR THE TIMEDIFF
      laglong[inext:inextnext,newstablecol]<-Timelags[,stablecol] #THIS IS FOR THE ANXIETY
      laglong[inext:inextnext,newvariablecol]<-Timelags[,imm] #THIS IS FOR THE ANXIETYLAG
      #if(!is.null(controlvariables)){
      #  for(cccccc in 1:length(controlvariables)){
      #    controldatacolumn=newvariablecol+cccccc
      #    laglong[inext:inextnext,controldatacolumn]<-Timelags[,imm]
      #  }
      #}
      if (!is.null(controllag)){
        for(cc in 1:length(controllag)){
          numberofcontrols=length(controllag)*i-length(controllag)
          lagcontrolcol=2*numberofvars+3+numberofcontrols
          #THIS NEXT PART IS THE PREVIOUS CODE USED TO CREATE NEW VARIABLES IN THE LONG FORMAT
          #  controllagmatrix<-matrix(NA,nrow=nrow(Timelags),ncol=(maxlags-controllag))
          #  Timelags<-data.frame(cbind(Timelags,controllagmatrix))
          #  laglong[inext:inextnext,lagcontrolcol]<-Timelags[,immm] #THIS IS THE CONTROL LAG CHANGE #ADDING A NEW VARIABLE #THIS MAY NOT WORK BECAUSE IT GOES BEYOND THE BOUNDS OF THE SCRIPT AND I'M NOT SURE IF IT WILL INSET NAs
          laglong[inext:inextnext,(lagcontrolcol+cc)]<-Timelags[,(stablecol+controllag[cc])]
        }
      }

    }
  }



  #MULTIVARIATE STACKING OF WIDE FORMAT (vector Timelags)
  numberofoutcomes=length(outcome) #COUNT NUMBER OF OUTCOMES

  #NOW STACK data FOR EACH OUTCOME
  blankdummyvars = matrix(NA, nrow = nrow(Timelags), ncol = numberofoutcomes)
  blankdummyvars = data.frame(blankdummyvars)
  names(blankdummyvars)[1:numberofoutcomes] = paste("Dummy", 1:numberofoutcomes, sep = "")
  Timelags = cbind(Timelags, blankdummyvars)
  namesofnewpredictorvariableswithnawide = NULL
  for(i in 1:numberofoutcomes){
    Timelagstemp = Timelags
    Timelagstemp[, paste("Dummy", i, sep = "")] = 1
    Timelagstemp[, paste("Dummy", as.vector(1:numberofoutcomes)[-i], sep = "")] = 0
    Timelagstemp$outcome = Timelagstemp[, paste(outcome[i], "lag0", sep = "")]
    for(ii in 1:length(differentialtimevaryingpredictors)){
      for(iii in 1:numberofoutcomes){
        temporarynewvarmat = matrix(NA, nrow = nrow(Timelags), ncol = numberofoutcomes)
        temporarynewvarmat = data.frame(temporarynewvarmat)
        #if(beforeblock){#NEW CODE JULY 10, 2016
        #  temporarynewvarmat=Timelagstemp[,paste("Dummy",iii,sep="")]*Timelagstemp[,paste(differentialtimevaryingpredictors[ii],"lag",(1:maxlags),sep="")] #OLD CODE
        #  names(temporarynewvarmat)[1:maxlags]=paste(differentialtimevaryingpredictors[ii],"lagon",outcome[iii],"lag",(1:maxlags),sep="") #OLD CODE
        #}else{#NEW CODE JULY 10, 2016
        temporarynewvarmat = Timelagstemp[, paste("Dummy", iii, sep = "")] * Timelagstemp[, paste(differentialtimevaryingpredictors[ii],
                                                                                                  "lag",
                                                                                                  (1:maxlags) * predictionsinterval,
                                                                                                  sep = "")] #July 6 2016, EDIT ADDED *time
        names(temporarynewvarmat)[1:maxlags] = paste(
          differentialtimevaryingpredictors[ii],
          "lagon",
          outcome[iii],
          "lag",
          (1:maxlags) * predictionsinterval,
          sep = ""
        ) #OLD CODE
        #}#NEW CODE JULY 10, 2016
        if (iii == 1) {
          temporarynewvarmatBACKUP = temporarynewvarmat
        } else{
          temporarynewvarmatBACKUP = cbind(temporarynewvarmatBACKUP, temporarynewvarmat)
        }
      }

      Timelagstemp = cbind(Timelagstemp, temporarynewvarmatBACKUP)
      namesofnewpredictorvariableswithnawide = unique(as.vector(c(
        namesofnewpredictorvariableswithnawide,
        paste(
          differentialtimevaryingpredictors[ii],
          "lagon",
          outcome[i],
          "lag",
          maxlags,
          sep = ""
        )
      )))
      namesofnewpredictorvariableswide = namesofnewpredictorvariableswithnawide[!is.na(namesofnewpredictorvariableswithnawide)]
    }
    if(i==1) {
      Timelagsdummy <- Timelagstemp
    } else{
      Timelagsdummy <- rbind(Timelagsdummy, Timelagstemp)
    }
  }

  ###END OF MULTIVARIATE STACKING FOR WIDE data


  laglongreduce <-
    laglong #na.exclude(laglong) TRYING TO KEEP THE NA data
  laglongreduce <- data.frame(laglongreduce)
  laglongreduce[, 3] <- as.numeric(laglongreduce[, 3])
  laglongreduce <- data.frame(laglongreduce)
  #NAMING COLUMNS
  names(laglongreduce)[1] <- "ID"
  names(laglongreduce)[2] <- "time"
  names(laglongreduce)[3] <- "timediff"
  for (i in 1:numberofvars){
    colnamestable = 2 + 2 * i
    colnamevar = 3 + 2 * i
    lagvarname = paste(input_list[i], "lag", sep = "")
    names(laglongreduce)[colnamestable] <- input_list[i]
    names(laglongreduce)[colnamevar] <- lagvarname
    if (!is.null(controllag)) {
      for (ii in 1:length(controllag)) {
        controllagname <-
          paste(input_list[i], "controllag", controllag[ii], sep = "")
        numberofcontrols = length(controllag) * i - length(controllag)
        controllagcolumn = 2 * numberofvars + 3 + numberofcontrols + ii
        names(laglongreduce)[controllagcolumn] <- controllagname
      }
    }
  }


  #2014.10.24
  #ADDING MULTIVARIATE OUTCOME
  numberofoutcomes = length(outcome) #COUNT NUMBER OF OUTCOMES
  #NOW STACK data FOR EACH OUTCOME
  blankdummyvars = matrix(NA, nrow = nrow(laglongreduce), ncol = numberofoutcomes)
  blankdummyvars = data.frame(blankdummyvars)
  names(blankdummyvars)[1:numberofoutcomes] = paste("Dummy", 1:numberofoutcomes, sep =
                                                      "")
  laglongreduce = cbind(laglongreduce, blankdummyvars)
  namesofnewpredictorvariableswithna = NA
  for (i in 1:numberofoutcomes) {
    laglongreducetemp = laglongreduce
    laglongreducetemp[, paste("Dummy", i, sep = "")] = 1
    laglongreducetemp[, paste("Dummy", as.vector(1:numberofoutcomes)[-i], sep =
                                "")] = 0
    laglongreducetemp$outcome = laglongreducetemp[, outcome[i]]
    for (ii in 1:length(differentialtimevaryingpredictors)) {
      temporarynewvarmat = matrix(NA, nrow = nrow(laglongreduce), ncol = numberofoutcomes)
      temporarynewvarmat = data.frame(temporarynewvarmat)
      temporarynewvarmat = laglongreducetemp[, paste("Dummy", as.vector(1:numberofoutcomes), sep =
                                                       "")] * laglongreducetemp[, paste(differentialtimevaryingpredictors[ii], "lag", sep =
                                                                                          "")]
      temporarynewvarmat = data.frame(temporarynewvarmat)
      names(temporarynewvarmat)[1:numberofoutcomes] = paste(differentialtimevaryingpredictors[ii],
                                                            "lagon",
                                                            outcome,
                                                            sep = "")
      laglongreducetemp = cbind(laglongreducetemp, temporarynewvarmat)
      namesofnewpredictorvariableswithna = unique(as.vector(c(
        namesofnewpredictorvariableswithna,
        paste(
          differentialtimevaryingpredictors[ii],
          "lagon",
          outcome[i],
          sep = ""
        )
      )))
      namesofnewpredictorvariables <-
        namesofnewpredictorvariableswithna[!is.na(namesofnewpredictorvariableswithna)]
    }
    if (i == 1) {
      laglongreducedummy = laglongreducetemp
    } else{
      laglongreducedummy = rbind(laglongreducedummy, laglongreducetemp)
    }
  }
  laglongreducedummybackup = laglongreducedummy
  laglongreducedummy <-
    laglongreducedummy[laglongreducedummy[, 3] > 0, ] #DELETING CONCURRENT TIMES
  #print(paste("predictionsend is equal to",predictionsend,"the number of rows in laglongreducedummy is ", nrow(laglongreducedummy)))
  if(!is.null(predictionsend)){
    laglongreducedummy <- laglongreducedummy[laglongreducedummy[,3]<=predictionsend,] #DELETING TIMES LESS THEN THE PREDICTIONS END
    #print(paste("the number of rows in laglongreducedummy is ", nrow(laglongreducedummy)))
  }
  laglongreducedummy<-laglongreducedummy
  #END OF data STACKING FOR MULTIVARIATE OUTCOMES

  datamanipulationoutlist=list("laglongreducedummy"=laglongreducedummy,"Timelagsdummy"=Timelagsdummy,"namesofnewpredictorvariables"=namesofnewpredictorvariables,"laglongmatrixlength"=laglongmatrixlength,"laglongreducedummybackup"=laglongreducedummybackup,"lengthcovariates"=lengthcovariates)

  #FIGURE OUT WHAT I WANT TO RETURN HERE!
  return(datamanipulationoutlist)
}
