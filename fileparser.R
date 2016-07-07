
con  <- file(file.choose(), open = "r")
alternateName=FALSE

#parse for alternative name scheme
oneLine <- readLines(con, n = 1, warn = FALSE)
trueFalse<- strsplit(oneLine,"=")[[1]]
if(trueFalse[[2]]=="TRUE")
{
  alternateName=TRUE
}

oneLine <- readLines(con, n=1, warn = FALSE)
print(oneLine)

#while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
#{
#  print(oneLine)
#} 

close(con)
