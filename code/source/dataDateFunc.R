
#date defining function
dataDateMkr <- function(dateVar=NA){
  #try get date from source if exist 
  args<-commandArgs(trailingOnly=TRUE) #pull from outside var when sourcing script
  #make globalDate if 
  dataDate<-if(length(args) > 0) {
    as.Date(args[1]) #if globalDate is NA & dateVar is not NA, set date in code with dateVar
  } else if(!is.na(dateVar)) {
    as.Date(dateVar)
  } else {
    Sys.Date()-1 #or sysDate -1
  }
  return(dataDate)
}