compositeField <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                      lty=1,col="black",lwd=1) {
  results <- composite.field(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                      lty=1,col="black",lwd=1)
  invisible(results)
}
