# Strips the strings by cutting off at the first space
# R.E. Benestad, Oslo, Norway, April 19.04.2002.
# met.no

strip<-function(in.str) {
  
  lfac<-FALSE                           # Set flag if we are dealing with a factor
                                        # object. Then the output is converted to
                                        # factor.
  if (is.factor(in.str)) { lfac <- TRUE }
  in.str<-as.character(in.str)
  
  strip <- in.str
  
# Go through list of a string array and remove the remainder of the string
# starting at the first space character...
  
  for (i in 1:length(strip)) {
    
# Leading spaces

    while (instring(" ",in.str[i])[1]==1) in.str[i] <- substr(in.str[i],2,nchar(in.str[i]))  

#  print(in.str)
    out.str<-in.str[i]

    c.str<-paste(unlist(strsplit(in.str[i],"")),sep="")
#    print(c.str)
    ispc<-pmatch(" ",c.str)
#    print(ispc)
    if (!is.na(ispc) & (ispc > 1)) {
      out.str<-substr(in.str[i],1,ispc-1)
    } 
  
    if (lfac) {
      out.str<-factor(out.str)
    }
    strip[i]<-out.str
  }
  strip
}
