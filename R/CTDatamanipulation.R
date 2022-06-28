datamanipulation = function(differentialtimevaryingpredictors = NULL,
                            outcome = NULL,
                            data = NULL,
                            ID = "ID",
                            Time = "time",
                            standardized = FALSE,
                            predictionsend = NULL) {


  # note - you probably don't want to do this!
  if(is.null(predictionsend)){
    warning("predictionsend using max time diff in dataset. Are you sure you want this?")
    predictionsend <- diff(data[nrow(data),Time], data[1,Time])
  }

  names_matrix = as.matrix(expand.grid(differentialtimevaryingpredictors, outcome))
  namesofnewpredictorvariables = c()
  for (i in 1:nrow(names_matrix)){
    addnames = paste(names_matrix[i,1],names_matrix[i,2],sep = "lagon" )
    namesofnewpredictorvariables = c(namesofnewpredictorvariables,addnames)
  }

  colnames = c(Time,"timediff","outcome",namesofnewpredictorvariables)
  datasep = matrix(ncol = length(colnames), nrow = 1)
  N = nrow(data)
  # if(predictionsend>N){
  #   predictionsend = N
  # }
  # for (i in 1:predictionsend) {
  #   for (j in 1:length(outcome)) {
  #     dataadd = matrix(ncol = length(colnames), nrow = (N-i))
  #     dataadd[,1] = data[,2][(i+1):N] # time of data need to start from 0
  #     dataadd[,2] = diff(data[,2],i)
  #     dataadd[,3] = data[,outcome[j]][(i+1):N]
  #
  #     for (k in 1:nrow(names_matrix)) {
  #       if(names_matrix[k,2] == outcome[j]){
  #         dataadd[,(k+3)] = data[,names_matrix[k,1]][1:(N-i)]
  #       } else{
  #         dataadd[,(k+3)] = 0
  #       }
  #
  #     }
  #     datasep = rbind(datasep,dataadd)
  #   }
  # }
  #
  for (i in 1:(N-1)) {
    time_diff = diff(data[,Time], lag = i)
    if(min(time_diff) <= predictionsend){
      for (j in 1:length(outcome)) {
        add_time_loc = which(time_diff <= predictionsend)  # find the location of previous time point which has the time diff less than the prediction end
        dataadd = matrix(ncol = length(colnames), nrow = length(add_time_loc)) # build the add-data matrix
        dataadd[,1] = data[,2][add_time_loc + i] # Fulfill the 1st column of add-data matrix, i.e. the time column
        dataadd[,2] = time_diff[add_time_loc] # Fulfill the 2nd column of add-data matrix, i.e. the time-diff column
        dataadd[,3] = data[,outcome[j]][add_time_loc + i] # Fulfill the 3rd column of add-data matrix, i.e. the outcome column
        for (k in 1:nrow(names_matrix)) {
          if(names_matrix[k,2] == outcome[j]){
            dataadd[,(k+3)] = data[,names_matrix[k,1]][add_time_loc] # Fulfill the value of Y#lagon Outcome, i.e., the previouse-i-lag value of this predictor
          } else{
            dataadd[,(k+3)] = 0
          }
        }
        datasep = rbind(datasep,dataadd)
      }
    }
  }
  laglongreducedummy = datasep[-1,] # remove the starting NA row
  colnames(laglongreducedummy) = colnames
  laglongreducedummy = as.data.frame(laglongreducedummy)
  lengthcovariates = length(namesofnewpredictorvariables) + 2
  return(list("laglongreducedummy" = laglongreducedummy, "namesofnewpredictorvariables" = namesofnewpredictorvariables,"lengthcovariates" = lengthcovariates ))
}




