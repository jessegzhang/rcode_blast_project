#3 test cases
#test<-"(((AID,BED),CAB),(DAD,EGG))"
#test<-"((ANC:.1,BED:.2):.15,(CAS:.1,DED:.05):.2)"
test<-"((ABE,BOG),((COD,DON),EGO))"

#phase out unecessary data in the newick notation
simplified<-gsub("[:.0-9]","",test)
myVars<-strsplit((gsub("[()]","",simplified)),"[,]")[[1]]

#initalizing list to store variables
varList <- list()
for (var in myVars)
{
  varList[[ var ]] <- 0
}

#splits the matrix into variables
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
    varList[[count]]<-priority
    count<-count+1
  }
}

#creates matrix 
finalMatrix <- matrix(0,nrow=length(varList),ncol=length(varList))
#names(varList) might be able to be subsitute by myVars 
colnames(finalMatrix)<-myVars
rownames(finalMatrix)<-myVars

#begin implementing calculation
#TODO: find a way to find pairs of variables
#3 cases
#2 variables with empty matricies = 1 variable with a matricie
#1 variable with empty 1 varibale with filled matrix = 1 variable 
# 2 variables with filled matrix


#generation of matrixList
#possible to fuse in previous for loop?
matrixList<-list()
for(i in 1:length(varList))
{
  matrixList[[i]]<- list(var=names(varList[i]), value=varList[[i]], varMatrix=matrix(0,nrow=length(varList),ncol=1))
  rownames(matrixList[[i]]$varMatrix)<-myVars
}

#finding highest priority and amount
varVal <- sapply(matrixList,"[[","value")
maxPriority<-varVal[which.max(abs(varVal))]


#ordering names in order
nameStore<-myVars[order(varVal, decreasing=TRUE)]
count<-0

#find appropriate indexes to modify
#TODO: maybe add a conditional statement to make sure numbers are even
for(i in 1:sum(varVal==maxPriority))
{
  print(match(nameStore[i],myVars))
}



