#3 test cases
#test<-"(((AID,BED),CAB),(DAD,EGG))"
#test<-"((ANC:.1,BED:.2):.15,(CAS:.1,DED:.05):.2)"
#test<-"(A,B)"
#test<-"((A,B),C)"
#test<-"((A,B),(C,D))"
#test<-"((ABE,BOG),((COD,DON),EGO))"
test<-"((A,B),((C,D),E))"
uniquevarID<-0
tempVar<-"tempVar"

#phase out unecessary data in the newick notation
simplified<-gsub("[:.0-9]","",test)
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
full_split<-gsub("[A-Z]"," ",simplified)
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
      for(i in myVars)
      {
        
        value<-matrixList[[indexList[[i]]]]$varMatrix[[i, 1]]
        if(value>0)
        {
          for(j in myVars)
          {
            secondValue<-matrixList[[indexList[[i+1]]]]$varMatrix[[i, 1]]
            if(secondValue>0)
            {
              finalMatrix[j,i]<-value+secondValue+2
              finalMatrix[i,j]<-value+secondValue+2
            }
          }
        }
      }
      
      #developing a new temp variable
      matrixList[[indexList[[i]]]]$value<-matrixList[[indexList[[i]]]]$value-1
      matrixList[[indexList[[i]]]]$varMatrix=matrixList[[indexList[[matrixVar]]]]$varMatrix
      
      matrixList[[indexList[[i]]]]$varMatrix<-matrixList[[indexList[[i]]]]$varMatrix+matrixList[[indexList[[i+1]]]]$varMatrix
      for(k in myVars)
      {
        value<-matrixList[[indexList[[i]]]]$varMatrix[[i, 1]]
        if(value>0)
        {
          finalMatrix[[k,1]]<-finalMatrix[[k,1]]+1
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
matrixCount<-0
#check first value
if(sum(matrixList[[1]]$varMatrix)>0)
{
  matrixCount<-matrixCount+1
  matrixVar<-1
  onlyVar<-matrixList[[2]]$var
}

#check second value
if(sum(matrixList[[2]]$varMatrix)>0)
{
  matrixCount<-matrixCount+1
  matrixVar<-2
  onlyVar<-matrixList[[1]]$var
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
        if(secondValue>0)
        {
          finalMatrix[j,i]<-value+secondValue+1
          finalMatrix[i,j]<-value+secondValue+1
        }
      }
    }
  }
}

#matrixList[[3]]<-NULL
#sum(matrixList[[3]]$varMatrix)
