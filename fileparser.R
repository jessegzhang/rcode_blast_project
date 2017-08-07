fileParser<-function(con)
{

  #TODO: add check for empty file
  
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
          varListKey[[1]]<-paste(keySplit[[2]],keySplit[[1]],sep="")
          keyList[[keyLength+1]]<-list(key=keySplit[[1]], varList=varListKey, length=1)
          keyLength<-keyLength+1
          #else we will be using the keyMatch index value and adding it to the varList for the key
          #and incrementing its length by 1
        } else {
          listLength<-keyList[[keyMatch]]$length
          keyList[[keyMatch]]$varList[[listLength+1]]<-paste(keySplit[[2]],keySplit[[1]],sep="")
          keyList[[keyMatch]]$length<-listLength+1
        }
        #replace a whitespace with a variable
        removeAltNamingScheme<-sub(" ",paste(keySplit[[2]],keySplit[[1]],sep=""),removeAltNamingScheme)
        print(removeAltNamingScheme)
      }
      #generates distance matrixes and input them into the distanceMatrixLIst
      distanceMatrixList[[DMLength+1]]<-grabDistanceMatrix(removeAltNamingScheme)
      DMLength<-DMLength+1
    }
    #if not using alternate naming scheme we generate distance matrixes normally
    else 
    {
      tempVarGrab<-simplifyNewick(oneLine)
      distanceMatrixList[[DMLength+1]]<-grabDistanceMatrix(tempVarGrab)
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
  
  #reduction of keyList
  tempKeyList<-vector(mode="list",length=keyLength)
  for(i in 1:keyLength)
  {
    tempVarListKey<-vector(mode="list", length=keyList[[i]]$length)
    for(j in 1:keyList[[i]]$length)
    {
      tempVarListKey[[j]]<-keyList[[i]]$varList[[j]]
    }
    tempKeyList[[i]]<-list(key=keyList[[i]]$key, varList=tempVarListKey)
  }
  keyList<-tempKeyList
  
  #print out for user confirmation
  #whether or not they like that list
  
  for(i in 1:keyLength)
  {
    cat("Key: ")
    cat(keyList[[i]]$key)
    cat("\nAssociated Taxa: ")
    for(j in keyList[[i]]$varList)
    {
      cat(j)
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
  #takes matrixReductionList and give it the associated entries denoted with 1
  for(i in 1:length(matrixReductionList))
  {
    matrixCount<-1
    for(j in colnames(matrixReductionList[[i]]$matrix))
    {
      taxaList <- sapply(keyList,"[[","varList")
      keyCount<-1
      for(k in taxaList)
      {
        keyMatch<-match(j, k)
        if(!is.na(keyMatch))
        {
          break
        }
        keyCount<-keyCount+1
      }
      matrixReductionList[[i]]$matrixKeyList[[keyCount]]$associatedMatrix[1,matrixCount]<-1
      matrixCount<-matrixCount+1
    }
  }
  
  
  
  #creates final phylogenetic matrix
  phyloMatrix <- matrix(0,nrow=length(keyList),ncol=length(keyList))
  colnames(phyloMatrix)<-sapply(keyList,"[[","key")
  rownames(phyloMatrix)<-sapply(keyList,"[[","key")
  #triple nested for loop
  for(i in seq.int(1,length(keyList)-1,1))
  {
    for(j in seq.int(i+1,length(keyList),1))
    {
      entryTotal<-0
      for(k in matrixReductionList)
      {
        #calculates entryCount and entrySum
        #takes first index and second index of the first two matrix in keylist and multiplies them
        entryCount<-sum(k$matrixKeyList[[i]]$associatedMatrix)*sum(k$matrixKeyList[[j]]$associatedMatrix)
        #t denotes the transpose of a matrix
        entrySum<-sum(k$matrixKeyList[[i]]$associatedMatrix%*%k$matrix%*%t(k$matrixKeyList[[j]]$associatedMatrix))
        
        if(entrySum==0)
        {
          stop("Matrix does not contain all species, exiting...")
        }
        entryTotal<-entryTotal+(entrySum/entryCount)
        #sums everything
      }
      #enter the entries in the phylogenetic matrix with the average of the averages
      #takes the average via length(matrixReductionlist)
      phyloMatrix[keyList[[i]]$key,keyList[[j]]$key]<-(entryTotal/length(matrixReductionList))
      phyloMatrix[keyList[[j]]$key,keyList[[i]]$key]<-(entryTotal/length(matrixReductionList))
    }
  }
  print(phyloMatrix)
  return(phyloMatrix)
}
#sets working directory for the other file
#setwd("C:/Users/jgzhang/Documents/Github/rcode_blast_project")
setwd("D:/GitHub/rcode_blast_project")
source("njst.R")
start.time <- Sys.time()

con  <- file(file.choose(), open = "r")

fileParser(con)

close(con)

end.time <- Sys.time()

time.taken <- end.time - start.time
