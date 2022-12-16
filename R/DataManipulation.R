datamanipulation = function(differentialtimevaryingpredictors = NULL,
                            outcome = NULL,
                            data = NULL,
                            ID = "ID",
                            Time = "time",
                            standardized = FALSE,
                            predictionsend = NULL) {

  # Find the time diff and the number of time points for each subject
  id_count = numeric(length(unique(data[,ID])))
  id_vec = unique(data[,ID])
  diffs = id_count
  numerpoints = id_count
  for (i in c(1:length(id_vec))){
    subset = data[which(data[,ID] == id_vec[i]),]
    diffs[i] = diff(range(subset[,Time]))
    numerpoints[i] = nrow(subset)
  }

  # Use the min diff to be the predictionsend if this argument is Null
  if(is.null(predictionsend)){
    warning("predictionsend using min time diff of subject within dataset. Are you sure you want this?")
    predictionsend = min(diffs)
  }
  # Create the name of covariates
  names_matrix = as.matrix(expand.grid(differentialtimevaryingpredictors, outcome))
  namesofnewpredictorvariables = c()

  for (i in 1:nrow(names_matrix)){
    addnames = paste(names_matrix[i,1],names_matrix[i,2],sep = "lagon" )
    namesofnewpredictorvariables = c(namesofnewpredictorvariables,addnames)
  }

  colnames = c(ID,Time,"timediff","outcome",namesofnewpredictorvariables)
  datasep = matrix(ncol = length(colnames), nrow = 1)

    for (ii in c(1:length(id_vec))){ # Stack data by ID
      for (i in c(1:(numerpoints[ii]-1))){ # Stack data by number of time points
    data_sub = data[which(data[,ID] == id_vec[ii]),]
    time_diff = diff(data_sub[,Time], lag = i)
    if(min(time_diff) <= predictionsend){
      for (j in 1:length(outcome)) {
        add_time_loc = which(time_diff <= predictionsend)  # find the location of previous time point which has the time diff less than the prediction end
        dataadd = matrix(ncol = length(colnames), nrow = length(add_time_loc)) # build the add-data matrix
        dataadd[,1] = id_vec[ii] # Fulfill the 1st column of add-data matrix, i.e. the ID column
        dataadd[,2] = data_sub[,Time][add_time_loc + i] # Fulfill the 2nd column of add-data matrix, i.e. the time column
        dataadd[,3] = time_diff[add_time_loc] # Fulfill the 3rd column of add-data matrix, i.e. the time-diff column
        dataadd[,4] = data_sub[,outcome[j]][add_time_loc + i] # Fulfill the 4th column of add-data matrix, i.e. the outcome column
        for (k in 1:nrow(names_matrix)) {
          if(names_matrix[k,2] == outcome[j]){
            dataadd[,(k+4)] = data_sub[,names_matrix[k,1]][add_time_loc] # Fulfill the value of Y#lagon Outcome, i.e., the previous-i-lag value of this predictor
          } else{
            dataadd[,(k+4)] = 0
          }
        }
        datasep = rbind(datasep,dataadd)
      }
    }
    }
  }

  laglongreducedummy = datasep[-1,] # remove the starting NA row
  colnames(laglongreducedummy) = colnames
  laglongreducedummy = as.data.frame(laglongreducedummy)
  lengthcovariates = length(namesofnewpredictorvariables) + 2
  return(list("laglongreducedummy" = laglongreducedummy,
              "namesofnewpredictorvariables" = namesofnewpredictorvariables,
              "lengthcovariates" = lengthcovariates ))
}

 test_all = datamanipulation(differentialtimevaryingpredictors = c("Y1","Y2"),outcome = c("Y1","Y2"),data = data, ID = "id", Time = "time",predictionsend = 20)


# data_test = rbind(data[1:3,], data[81:84,])
# test = datamanipulation(differentialtimevaryingpredictors = c("Y1","Y2"),outcome = c("Y1","Y2"),data = data_test, ID = "id", Time = "time",predictionsend = 20)
# test_test = datamanipulation_test(differentialtimevaryingpredictors = c("Y1","Y2"),outcome = c("Y1","Y2"),data = data_test[1:3,], ID = "id", Time = "time",predictionsend = 20)
