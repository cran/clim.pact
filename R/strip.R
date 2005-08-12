# Strips the strings by cutting off at the first space
# R.E. Benestad, Oslo, Norway, April 19.04.2002.
# met.no

strip<-function(in.str) {
  
  lfac<-FALSE                           # Set flag if we are dealing with a factor
                                        # object. Then the output is converted to
                                        # factor.
  if (is.factor(in.str)) { lfac <- TRUE }
  in.str<-as.character(in.str)
  
# Leading spaces
  
  while (instring(" ",in.str)[1]==1) in.str <- substr(in.str,2,nchar(in.str))  

#  print(in.str)
  out.str<-in.str

# Go through list of a string array and remove the remainder of the string
# starting at the first space character...
  
  for (is in 1:length(in.str)) {
    c.str<-paste(unlist(strsplit(in.str[is],"")),sep="")
#    print(c.str)
    ispc<-pmatch(" ",c.str)
#    print(ispc)
    if (!is.na(ispc) & (ispc > 1)) {
      out.str[is]<-substr(in.str[is],1,ispc-1)
    } 
  }
  if (lfac) {
    out.str<-factor(out.str)
  }
  strip<-out.str
  strip
}
