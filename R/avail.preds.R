avail.preds <- function(direc="data") {
  dir.0 <- getwd()
  if (file.exists(direc)) setwd(direc)
  avail.preds <- list.files(pattern=".Rdata")
  avail.preds <- avail.preds[grep("eof",avail.preds)]
  setwd(dir.0)
  avail.preds
}
