avail.ds <- function(direc="output") {
  dir.0 <- getwd()
  if (file.exists(direc)) setwd(direc)
  avail.ds <- list.files(pattern=".Rdata")
  avail.ds <- avail.ds[grep("ds",avail.ds)]
  setwd(dir.0)
  avail.ds
}
