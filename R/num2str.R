# Converts a number to string and formats with requested
# decimal points
# R.E. Benestad

num2str <- function(x,dec=2,f.width=NULL,d.point=".") {
  num2str <- rep(" ",length(x))
  for (i in 1:length(x)) {
    x.r <- as.character(round(x[i],dec))
    dot <- instring(".",x.r)
    if (dot==0) {
      x.r <- paste(x.r,".",sep="")
      dot <- instring(".",x.r)
    }
    zeros <- ""
#    print(x.r)
#    print(c(nchar(x.r),dot,NA,nchar(x.r)-dot+1,dec))
    if (nchar(x.r)-dot+1 <= dec) {
      for (ii in (nchar(x.r)-dot+1):dec) {
        zeros <- paste(zeros,"0",sep="")
      }
    }
    y <- paste(x.r,zeros,sep="")
    num2str[i] <- y
  }
  num2str
}
