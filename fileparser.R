#needs to be changed appropriately
setwd("~/GitHub/rcode_blast_project")
source("njst.R")

con  <- file(file.choose(), open = "r")

alternateName=FALSE

keyList<-vector(mode="list", lengt=1000)
keyLength<-1

#parse for alternative name scheme
oneLine <- readLines(con, n = 1, warn = FALSE)
trueFalse<- strsplit(oneLine,"=")[[1]]
if(trueFalse[[2]]=="TRUE")
{
  alternateName=TRUE
}

#throw out the second //
oneLine <- readLines(con, n=1, warn = FALSE)

distanceMatrixList <- vector(mode = "list", length = 1000)

DMLength<-1
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0 && oneLine!="//")
{
  if(alternateName==TRUE)
  {
    tempVarGrab<-simplifyNewick(oneLine)
    matrixVars<-strsplit((gsub("[()]","",tempVarGrab)),"[,]")[[1]]
    for(i in matrixVars)
    {
      keySplit<-strsplit(i,"_")[[1]]
      print(keySplit[[1]])
      for(j in 1:DMLength)
      {
        
      }
    }
  }
  else 
  {
  distanceMatrixList[[DMLength]]<-grabDistanceMatrix(oneLine)
  DMLength<-DMLength+1
  }
} 



if(alternateName==FALSE)
{
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)))
  {
    #collapses whitespace
    splitString<-gsub(" ","",oneLine)
    splitString <- strsplit(splitString, ":")[[1]]
    varListKey<-strsplit(splitString[[2]], ",")[[1]]
    keyList[[keyLength]]<-list(key=splitString[[1]], varList=varListKey, varLength=length(varListkey))
    keyLength<-keyLength+1
  }
}

#add associated print menu here
close(con)
