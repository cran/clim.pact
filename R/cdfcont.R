# Obtains the name of variables in a netcdf File
#
# R.E. Benestad, 23.09.2003

cdfcont <- function(filename) {
  system(paste("ncdump -h  ",filename," > cdfcont.txt",sep=""),intern=T)
  cdfhead <- readLines("cdfcont.txt")
  cdfvars <- cdfhead[c(grep("float",cdfhead),
                       grep("short",cdfhead),
                       grep("double",cdfhead))]
  cdfdims <- cdfvars
  for (i in 1:length(cdfvars)) {
    i1 <- instring(" ",cdfvars[i])
    i2 <- instring("(",cdfvars[i])
    i3 <- instring(")",cdfvars[i])
    cdfdims[i] <- substr(cdfvars[i],i2[1]+1,i3[1]-1)
    cdfvars[i] <- substr(cdfvars[i],i1[1]+1,i2[1]-1)
  }
#  print(cdfvars)
#  print(cdfdims)
  system("rm -f cdfcont.txt",intern=T)
  content <- list(vars=cdfvars,dims=cdfdims)
  invisible(content)
}
