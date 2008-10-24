# R.E. Benestad, met.no, Oslo, Norway 09.10.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------



newFig <- function() {
#   sinfo <- Sys.info()
#   if (lower.case(sinfo[1])=="linux") x11() else
#                                      windows()
#   if (class(options()$device)=="character") {
#     device.type <- options()$device[1]
#     dev <- paste(device.type,"()",sep="") 
#     #print(paste("newFig: options()$device=",dev))
#     if ((dev!="none()") & (dev!="bitmap()")) eval(parse(text=dev)) else
#     if (dev=="bitmap()") {
#       while (dev.cur() > 1) dev.off()
#       bitmap(file="newFig.jpg",type="jpeg")
#     }
#   } else {
#     print(paste("You seem to run R on",sinfo[1],"which doesn't run smoothly with clim.pact -",
#                 "maybe you should consider using Linux?"))
#   }
  dev.new()
}
