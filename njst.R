#3 test cases
#test<-"(((AID,BED),CAB),(DAD,EGG))"
#test<-"((AN3C:.1,B1ED:.2)123:.15,(CAS2:.1,DE4D:.05)125:.2)"
#test<-"(A,B)"
#test<-"((A,B),C)"
#test<-"((A,B),(C,D))"
#test<-"((ABE,BOG),((COD,DON),EGO))"
#test<-"((A,B),((C,D),E))"
#test<-"(((A,B),(C,D)),(E,F))"
#test<-"((A,B),C,D)"
#test<-"((A,B),(C,D),E)"
#test<-"((A,B),(C,D),(E,F))"
test<-"(((A,B),(C,D)),(E,F),(G,H))"
uniquevarID<-0
tempVar<-"tempVar"

#phase out unecessary data in the newick notation
test_split<- strsplit(test,"")[[1]]

count<-1
numCheck<-FALSE
for (char in test_split)
{
  
  if(char==":")
  {
    numCheck<-TRUE
  }
  
  if(char==","||char==")"||char=="(")
  {
    numCheck<-FALSE  
  }
  
  if(numCheck==TRUE)
  {
    test_split <- test_split[-count]
  }
  else
  {
    count<-count+1 
  }
}

#used to remove extra numbers from t he original pattern recognizer, may be combined in the future if possible
count<-1
numCheck<-FALSE
for (char in test_split)
{

  if(char!="0"&&char!="1"&&char!="2"&&char!="3"&&char!="4"&&char!="5"&&char!="6"&&char!="7"&&char!="8"&&char!="9")
  {
    numCheck<-FALSE
  }
  if(char==")")
  {
    numCheck<-TRUE
  }
  
  if(count>length(test_split))
  {
    break
  }
  if(numCheck==TRUE&&test_split[count]!=","&&test_split[count]!=")")
  {
    test_split<-test_split[-count]
  }
  else
  {
    count<-count+1
  }
}


simplified<-paste(test_split, collapse="")

myVars<-strsplit((gsub("[()]","",simplified)),"[,]")[[1]]

#initalizing list to store variables
#TODO: Preallocate size of list
matrixList <- vector(mode="list",length=length(myVars))
count<-1
for (vars in myVars)
{
  matrixList[[count]]<- list(var=vars, value=0, varMatrix=matrix(0,nrow=length(myVars),ncol=1))
  rownames(matrixList[[count]]$varMatrix)<-myVars
  count<-count+1
}

#splits the newick notation into variables
full_split<-gsub("[A-Z0-9]"," ",simplified)
full_split<-gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", full_split, perl=TRUE)
full_split <- strsplit(full_split, "")[[1]]

#establishes variables
priority<-0
count<-1
#places priority on variables
for(i in full_split)
{
  if(i=="(")
  {
    priority<-priority+1
  }
  if(i==")")
  {
    priority<-priority-1
  }
  if(i==" ")
  {
    matrixList[[count]]$value<-priority
    count<-count+1
  }
}

#creates matrix 
finalMatrix <- matrix(0,nrow=length(myVars),ncol=length(myVars))
#names(varList) might be able to be subsitute by myVars 
colnames(finalMatrix)<-myVars
rownames(finalMatrix)<-myVars

while(length(matrixList)>2)
{
  #begin implementing calculation
  #TODO: find a way to find pairs of variables
  #3 cases
  #2 variables with empty matricies = 1 variable with a matricie
  #1 variable with empty 1 varibale with filled matrix = 1 variable 
  # 2 variables with filled matrix
  
  
  #finding highest priority and amount
  varVal <- sapply(matrixList,"[[","value")
  maxPriority<-varVal[which.max(abs(varVal))]
  
  
  #ordering names in order
  nameStore<-myVars[order(varVal, decreasing=TRUE)]
  count<-0
  
  


  
  #find appropriate indexes to modify
  #TODO: maybe add a conditional statement to make sure numbers are even
  #list preallocated to be more efficient
  indexList<-vector("list",100)
  for(i in 1:sum(varVal==maxPriority))
  {
    indexList[[count+1]]<-match(nameStore[i],myVars)
    count<-count+1
  }
  
  #conditional if statement htat breaks out of a while loop if
  #matrixList has 3 variables and they all have the same remaining priority
  if(length(matrixList)==3&&count==3)
  {
    break
  }
  
  #handles two variables at a time
  #conditional may need to be added for odd ones
  #three cases in here
  for(i in seq.int(1,count,2))
  { 
    matrixCount<-0
    #check first value
    if(sum(matrixList[[indexList[[i]]]]$varMatrix)>0)
    {
      matrixCount<-matrixCount+1
    }
    
    #check second value
    if(sum(matrixList[[indexList[[i+1]]]]$varMatrix)>0)
    {
      matrixCount<-matrixCount+1
    }
    #establish a tempvariable
    newVar<- c(tempVar,uniquevarID)
    newVar<-paste(newVar, collapse="")
    uniquevarID<-uniquevarID+1
    
    #store variable names
    firstVar<-matrixList[[indexList[[i]]]]$var
    secondVar<-matrixList[[indexList[[i+1]]]]$var
    #three cases below
    #two variables
    if(matrixCount==0)
    {
      
      
      #place variables on the finalMatrix
      finalMatrix[firstVar,secondVar]<-2
      finalMatrix[secondVar,firstVar]<-2
      
      #modifying value to be a temp variable
      matrixList[[indexList[[i]]]]$value<-matrixList[[indexList[[i]]]]$value-1
      matrixList[[indexList[[i]]]]$varMatrix[firstVar,1]<-1
      matrixList[[indexList[[i]]]]$varMatrix[secondVar,1]<-1
      matrixList[[indexList[[i]]]]$var<-newVar
      
      #delete the entry that wasn't the temp variable
      #needs to be moved outside of for loop
    }
    
    #one variable, one matrix
    if(matrixCount==1)
    {
      #establish the proper variables
      #check for which variable is the matrix and which is not
      if(sum(matrixList[[indexList[[i]]]]$varMatrix)>0)
      {
        matrixVar<-i
        onlyVar<-secondVar
      } 
      else
      {
        matrixVar<-i+1
        onlyVar<-firstVar
      }
      
      
      for(j in myVars)
      {
        value<-matrixList[[indexList[[matrixVar]]]]$varMatrix[[j, 1]]
        if(value>0)
        {
          finalMatrix[onlyVar,j]<-value+2
          finalMatrix[j,onlyVar]<-value+2
        }
      }
      
      #modifying new temp variable value
      matrixList[[indexList[[i]]]]$value<-matrixList[[indexList[[i]]]]$value-1
      matrixList[[indexList[[i]]]]$varMatrix=matrixList[[indexList[[matrixVar]]]]$varMatrix
      
      #developing new matrix
      for(k in myVars)
      {
        value<-matrixList[[indexList[[i]]]]$varMatrix[k,1]
        if(value>0)
        {
          matrixList[[indexList[[i]]]]$varMatrix[k,1]<-matrixList[[indexList[[i]]]]$varMatrix[k,1]+1
        }
      }
      matrixList[[indexList[[i]]]]$varMatrix[onlyVar,1]<-1
      matrixList[[indexList[[i]]]]$var<-newVar
    }
    
    #two matrix variables
    #most complex computation
    if(matrixCount==2)
    {
      #fill in matrix values
      for(ii in myVars)
      {
        
        value<-matrixList[[indexList[[i]]]]$varMatrix[[ii, 1]]
        if(value>0)
        {
          for(jj in myVars)
          {
            secondValue<-matrixList[[indexList[[i+1]]]]$varMatrix[[jj, 1]]
            if(secondValue>0 && ii!=jj)
            {
              finalMatrix[jj,ii]<-value+secondValue+2
              finalMatrix[ii,jj]<-value+secondValue+2
            }
          }
        }
      }
      
      #developing a new temp variable
      matrixList[[indexList[[i]]]]$value<-matrixList[[indexList[[i]]]]$value-1
      matrixList[[indexList[[i]]]]$varMatrix<-matrixList[[indexList[[i]]]]$varMatrix+matrixList[[indexList[[i+1]]]]$varMatrix
      
      for(k in myVars)
      {
        value<-matrixList[[indexList[[i]]]]$varMatrix[[k, 1]]
        if(value>0)
        {
          matrixList[[indexList[[i]]]]$varMatrix[[k,1]]<- matrixList[[indexList[[i]]]]$varMatrix[[k,1]]+1
        }
      }
    }
    
  } 
  
  #deletion of redundant variables after phase
  for(i in seq.int(count,1,-2))
  {
    matrixList[[indexList[[i]]]]<-NULL
  }
}  

#Last 2 variables

if(length(matrixList)==2)
{
  matrixCount<-0
  #check first value
  if(sum(matrixList[[1]]$varMatrix)>0)
  {
    matrixCount<-matrixCount+1
    onlyVar<-matrixList[[2]]$var
    matrixVar<-1
  }
  
  #check second value
  if(sum(matrixList[[2]]$varMatrix)>0)
  {
    matrixCount<-matrixCount+1
    onlyVar<-matrixList[[1]]$var
    matrixVar<-2
  }
  
  firstVar<-matrixList[[1]]$var
  secondVar<-matrixList[[2]]$var
  
  #final case for 2 variable
  if(matrixCount==0)
  {
    finalMatrix[firstVar,secondVar]<-1
    finalMatrix[secondVar,firstVar]<-1
  }
  
  #final case for 1 variable and 1 matrix
  if(matrixCount==1)
  {
    for(i in myVars)
    {
      
      value<-matrixList[[matrixVar]]$varMatrix[[i, 1]]
      if(value>0)
      {
        finalMatrix[onlyVar,i]<-value+1
        finalMatrix[i,onlyVar]<-value+1
      }
    }
  }
  
  #final case for 2 matrixs
  if(matrixCount==2)
  {
    for(i in myVars)
    {
      
      value<-matrixList[[1]]$varMatrix[[i, 1]]
      if(value>0)
      {
        for(j in myVars)
        {
          secondValue<-matrixList[[2]]$varMatrix[[j,1]]
          if(secondValue>0 && i!=j)
          {
            finalMatrix[j,i]<-value+secondValue+1
            finalMatrix[i,j]<-value+secondValue+1
          }
        }
      }
    }
  }
}


if(length(matrixList)==3)
{
  matrixIndexes<-vector("list",3)
  
  #counts matrix and stores associated amount of matrixs
  matrixCount<-0
  for(i in 1:length(matrixList))
  {
    if(sum(matrixList[[i]]$varMatrix>0))
    {
      matrixCount<-matrixCount+1
      matrixIndexes[[matrixCount]]<-i
    }
  }
  

  #case for 0 matrix in 3 remaining taxa
  if(matrixCount==0)
  {
    firstVar<-matrixList[[1]]$var
    secondVar<-matrixList[[2]]$var
    thirdVar<-matrixList[[3]]$var
    
    finalMatrix[firstVar,secondVar]<-1
    finalMatrix[secondVar,firstVar]<-1
    finalMatrix[firstVar, thirdVar]<-1
    finalMatrix[thirdVar,firstVar]<-1
    finalMatrix[secondVar, thirdVar]<-1
    finalMatrix[thirdVar, secondVar]<-1
  }
  
  #case for 1 matrix in 3 remaining taxa
  if(matrixCount==1)
  {
    oneMatrix<-6-matrixIndexes[[1]]
    matrixVar<-matrixIndexes[[1]]
    #if statement to associate the proper non matrix variables
    if(oneMatrix==5)
    {
      firstVar<-matrixList[[2]]$var
      secondVar<-matrixList[[3]]$var
    }
    if(oneMatrix==4)
    {
      firstVar<-matrixList[[1]]$var
      secondVar<-matrixList[[3]]$var
    }
    if(oneMatrix==3)
    {
      firstVar<-matrixList[[1]]$var
      secondVar<-matrixList[[2]]$var
    }
    for(i in myVars)
    {
      value<-matrixList[[matrixVar]]$varMatrix[[i, 1]]
      if(value>0)
      {
        finalMatrix[firstVar,i]<-value+2
        finalMatrix[i,firstVar]<-value+2
        finalMatrix[secondVar,i]<-value+2
        finalMatrix[i,secondVar]<-value+2
      }
    }
    finalMatrix[firstVar,secondVar]<-2
    finalMatrix[secondVar,firstVar]<-2
  }
  
  #Case for 2 matrix in 3 remaining taxa
  if(matrixCount==2)
  {

    #Determines the non matrix by taking note of the indicies 
    oneMatrix<-6-matrixIndexes[[1]]-matrixIndexes[[2]]
    firstVar<-matrixList[[oneMatrix]]$var
    firstMatrix<-matrixIndexes[[1]]
    secondMatrix<-matrixIndexes[[2]]
    #combine two matrix temporarily to calculate distance to the one variable
    matrixCombined<-matrixList[[firstMatrix]]$varMatrix+matrixList[[secondMatrix]]$varMatrix
    for(i in myVars)
    {
      
      value<-matrixList[[firstMatrix]]$varMatrix[[i, 1]]
      if(value>0)
      {
        for(j in myVars)
        {
          secondValue<-matrixList[[secondMatrix]]$varMatrix[[j,1]]
          if(secondValue>0 && i!=j)
          {
            finalMatrix[j,i]<-value+secondValue+2
            finalMatrix[i,j]<-value+secondValue+2
          }
        }
      }
      
      value<-matrixCombined[[i,1]]
      if(value>0)
      {
        finalMatrix[i,firstVar]<-value+2
        finalMatrix[firstVar,i]<-value+2
      }
    }
    
  }
  
  #Case for 3 matrix in 3 remaining taxa
  if(matrixCount==3)
  {
    #combine two matrix to add for calculation in the third
    matrixCombined<-matrixList[[1]]$varMatrix+matrixList[[2]]$varMatrix
    
    for(i in myVars)
    {
      #calculate between two matrix
      value<-matrixList[[1]]$varMatrix[[i, 1]]
      if(value>0)
      {
        for(j in myVars)
        {
          secondValue<-matrixList[[2]]$varMatrix[[j,1]]
          if(secondValue>0 && i!=j)
          {
            finalMatrix[j,i]<-value+secondValue+2
            finalMatrix[i,j]<-value+secondValue+2
          }
        }
      }
      
      #calculate between the last two matrix
      value<-matrixCombined[[i,1]]
      if(value>0)
      {
        for(j in myVars)
        {
          secondValue<-matrixList[[3]]$varMatrix[[j,1]]
          if(secondValue>0 && i!=j)
          {
            finalMatrix[j,i]<-value+secondValue+2
            finalMatrix[i,j]<-value+secondValue+2
          }
        }
      }
    }
  }
}

