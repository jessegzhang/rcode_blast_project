#2 test cases
test<-"(((AID,BED),CAB),(DAD,EGG))"
#test<-"((ANC:.1,BED:.2):.15,(CAS:.1,DED:.05):.2)"

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
matrixtest <- matrix(0,nrow=length(varList),ncol=length(varList))
colnames(matrixtest)<-names(varList)
rownames(matrixtest)<-names(varList)

#begin implementing calculation
#TODO: find a way to find pairs of variables
#3 cases

https://stackoverflow.com/questions/28056727/obtain-key-from-value-in-r
