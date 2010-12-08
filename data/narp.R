as.matrix(read.table(file="narp.dat")) -> narp
d.narp <- dim(narp)
narp <- as.numeric(narp)
dim(narp) <- d.narp
rm(d.narp)

