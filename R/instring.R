instring <- function(c,target,case.match=TRUE) {
  l <- nchar(target)
  if (!case.match) {
    c <- lower.case(c)
    target <- lower.case(target)
  }
  pos <- 0
  for (i in 1:l) {
    tst <- substr(target,i,i)
   if (tst==c) pos <- c(pos,i)
   }
  if (length(pos) > 1) pos <- pos[-1]
  pos
}
