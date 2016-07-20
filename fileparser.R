#needs to be changed appropriately
setwd("~/GitHub/rcode_blast_project")
source("njst.R")

con  <- file(file.choose(), open = "r")

#perhaps check for empty file?

alternateName=FALSE

keyList<-vector(mode="list", length=1000)
keyLength<-0

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
#sets initial distance matrix length
DMLength<-0

#iterates through each of the lines and handles each line differently if alternate naming scheme is true
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0 && oneLine!="//")
{
  #alternate naming scheme
  if(alternateName==TRUE)
  {
    tempVarGrab<-simplifyNewick(oneLine)
    
    #grab matrixVars for the incoming for loop
    matrixVars<-strsplit((gsub("[()]","",tempVarGrab)),"[,]")[[1]]
    
    #remove each var into whitespace to be replaced by simplified var
    removeAltNamingScheme<-gsub("[A-Z0-9a-z]"," ",tempVarGrab)
    removeAltNamingScheme<-gsub("_"," ",removeAltNamingScheme)
    removeAltNamingScheme<-gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", removeAltNamingScheme, perl=TRUE)
    
    #iterates through each variable, looks at the value before the underscore and sees
    #if the entry is in keyList if not, generates an entry
    for(i in matrixVars)
    {
      #splits so keySplit is split into two indexes via the underscore
      keySplit<-strsplit(i,"_")[[1]]
      
      #grabs all key values
      varVal <- sapply(keyList,"[[","key")
      
      #attempts to match the first index(key) with the keyList
      #becomes NA if the key is not found
      keyMatch<-match(keySplit[[1]], varVal)
      
      #if the key is not found then it generates an entry
      if(is.na(keyMatch))
      {
        varListKey<-vector(mode="list", length=100)
        varListKey[[1]]<-keySplit[[2]]
        keyList[[keyLength+1]]<-list(key=keySplit[[1]], varList=varListKey, length=1)
        keyLength<-keyLength+1
        #else we will be using the keyMatch index value and adding it to the varList for the key
        #and incrementing its length by 1
      } else {
        listLength<-keyList[[keyMatch]]$length
        keyList[[keyMatch]]$varList[[listLength+1]]<-keySplit[[2]]
        keyList[[keyMatch]]$length<-listLength+1
      }
      #replace a whitespace with a variable
      removeAltNamingScheme<-sub(" ",keySplit[[2]],removeAltNamingScheme)
    }
    #generates distance matrixes and input them into the distanceMatrixLIst
    distanceMatrixList[[DMLength+1]]<-grabDistanceMatrix(removeAltNamingScheme)
    DMLength<-DMLength+1
  }
  #if not using alternate naming scheme we generate distance matrixes normally
  else 
  {
  distanceMatrixList[[DMLength+1]]<-grabDistanceMatrix(oneLine)
  DMLength<-DMLength+1
  }
} 

#at the end if the alternate name scheme is not being used we will have to make a new keyList based on
#supplied data
#e.g. H: AID, BED
if(alternateName==FALSE)
{
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)))
  {
    #collapses whitespace
    splitString<-gsub(" ","",oneLine)
    splitString <- strsplit(splitString, ":")[[1]]
    varListKey<-strsplit(splitString[[2]], ",")[[1]]
    listLength<-length(varListKey)
    keyList[[keyLength+1]]<-list(key=splitString[[1]], varList=varListKey, length=listLength)
    keyLength<-keyLength+1
  }
}

#print out for user confirmation
#whether or not they like that list

for(i in 1:keyLength)
{
  cat("Key: ")
  cat(keyList[[i]]$key)
  cat("\nAssociated Variables: ")
  for(j in 1:keyList[[i]]$length)
  {
    cat(keyList[[i]]$varList[[j]])
    cat(" ")
  }
  cat("\n\n")
}

matrixReductionList<-vector(mode="list", length=DMLength)
#establish list for reduction
for(i in 1:DMLength)
{
  matrixVarKey<-vector(mode="list", length=keyLength)
  tempMatrix=matrix(0, nrow=1, ncol=nrow(distanceMatrixList[[i]]))
  for(j in 1:keyLength)
  {
    matrixVarKey[[j]]<-list(key=keyList[[j]]$key, associatedMatrix=tempMatrix)
  }
  matrixReductionList[[i]]<-list(matrix=distanceMatrixList[[i]], matrixKeyList=matrixVarKey)
}

#calculation of the matrix here
for(i in matrixReductionList)
{
  print(i$matrix)
}

#add associated print menu here
close(con)