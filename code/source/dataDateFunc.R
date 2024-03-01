
#date defining function
dataDateMkr <- function(dateVar=NA){
  #try get date from source if exist 
  globalDate<-commandArgs(trailingOnly=TRUE)[1] #pull from outside var when sourcing script
  #make globalDate if 
  globalDate<-if(is.na(globalDate) & !is.na(dateVar)){
            as.Date(dateVar) #if globalDate is NA & dateVar is not NA, set date in code with dateVar
            }else{
            if(exists("globalDate")){
              as.Date(globalDate)
              }else{NA}
            } 
  dataDate<-if(!is.na(globalDate)){
                globalDate-1 #make dataDate from globalDate
                }else{
                Sys.Date()-1 #or sysDate -1
                }  
  return(dataDate)
  }