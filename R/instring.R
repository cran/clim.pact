instring <- function(c,target,case.match=TRUE) {
  l <- nchar(target)
  if (!case.match) {
    c <- lower.case(c)
    target <- lower.case(target)
  }

  nc <- nchar(c)
  
  if (nc ==1) {
# c is one character: use the old routine
    pos <- 0
    for (i in 1:l) {
      tst <- substr(target,i,i)
     if (tst==c) pos <- c(pos,i)
     }
    if (length(pos) > 1) pos <- pos[-1]
  } else { 
# c is a string: use new routine
#print(c)
    spos <- rep(NA,nc*l); dim(spos) <- c(nc,l)
    for (j in 1:nc) {
      a <- instring(substr(c,j,j),target)
      if (length(a) > 0) {
        if (j > 1) {
#print(paste("a=",a))
#print(paste("spos[j-1,]=",spos[j-1,]))
          a.match <- is.element(a,spos[j-1,]+1)
          a <- a[a.match]
          p.match <- is.element(spos[j-1,],a-1)
#print(paste("a.match:",sum(a.match),"  p.match:",sum(p.match)))
          spos <- spos[,p.match]
          dim(spos) <- c(nc,sum(p.match))
        }
#print(paste("a=",a))
#print(dim(spos))
        if (length(a) < dim(spos)[2]) {
          spos <- spos[,1:length(a)]
          dim(spos) <- c(nc,length(a))
#print(dim(spos))
        }
        spos[j,] <- a
      } else spos[j,] <- 0
    }
    ns <- dim(spos)[2]
#print(paste("spos=",spos,"  ns=",ns))
    if (ns > 0) {
      pos <- rep(0,ns)
      for (i in 1:ns) {
        b <- c(diff(spos[,i]),1)
        i1 <- is.element(b,1)
        b[!i1] <- 0
        if (sum(b) == nc) pos[i] <- spos[1,i]
      }
    } else {pos <- 0; pos <- pos[-1]}
  } 
  pos  
}
