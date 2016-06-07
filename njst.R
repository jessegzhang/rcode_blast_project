#3 test cases
#test<-"(((AID,BED),CAB),(DAD,EGG))"
#test<-"((ANC:.1,BED:.2):.15,(CAS:.1,DED:.05):.2)"
test<-"((ABE,BOG),((COD,DON),EGO))"
uniquevarID<-0

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
for(i in seq.int(1,count,2)) { 
  
} 
matrixList[[3]]<-NULL
#sum(matrixList[[3]]$varMatrix)

