compositeField <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                      lty=1,col="black",lwd=1,main=NULL,sub=NULL) {
  results <- composite.field(x,y,lsig.mask=lsig.mask,sig.lev=sig.lev,s=s,mon=mon,
                      lty=lty,col=col,lwd=lwd,main=main,sub=sub)
  invisible(results)
}
