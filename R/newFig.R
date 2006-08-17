# R.E. Benestad, met.no, Oslo, Norway 09.10.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------



newFig <- function() {
   dev <- paste(options()$device,"()",sep="")
   #print(paste("newFig: options()$device=",dev))
   if ((dev!="none()") & (dev!="bitmap()")) eval(parse(text=dev)) else
   if (dev=="bitmap()") {
   while (dev.cur() > 1) dev.off()
   bitmap(file="newFig.jpg",type="jpeg")
 }
 }
