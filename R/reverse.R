# Produces a vector with the reversed order to that of a sort
# call
# R.E. Benestad
reverse.sort <- function(x) {
  reverse <- -1*sort(-1*x)
  reverse
}

reverse <- function(x) {
  reverse  <- x
  for (i in 1:length(x)) reverse[i]  <-  x[length(x)-i+1]
  reverse
}
