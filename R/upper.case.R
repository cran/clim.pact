# Returns the upper case of a string;
# Test:
# > upper.case("abcdefghijklemnoprstuvwxyz1234567890")
# [1] "ABCDEFGHIJKLEMNOPRSTUVWXYZ1234567890"
# > upper.case("ABCDEFGHIJKLEMNOPRSTUVWXYZ1234567890")
# [1] "ABCDEFGHIJKLEMNOPRSTUVWXYZ1234567890"
# R.E. Benestad (REB), DNMI (met.no), 08.01.2002
#                REB:  03.05.2002 - modified to hande string arrays

upper.case <- function(u.case) {

  lfac<-FALSE                           # Set flag if we are dealing with a factor
                                        # object. Then the output is converted to
                                        # factor.
  
  if (is.null(u.case)) return()
  if (is.factor(u.case)) { lfac <- TRUE }
  if ( (!is.character(u.case)) & (!is.factor(u.case)) ) return()  
  str<-as.character(u.case)
#  print(min(nchar(str)))
  if ( (min(nchar(str))==0) & (is.null(str)) ) return()
#  print(str)  
  upper.case<-str

  for (is in 1:length(str)) {
    nc<-nchar(str[is])
    upper.case[is]<-""
    for (ic in 1:nc) {
      sstr<-substr(str[is],ic,ic)
      u.l<-switch(as.character(sstr),
                      a="A",b="B",c="C",d="D",f="F",e="E",g="G",h="H",i="I",
                      j="J",k="K",l="L",m="M",n="N",o="O",p="P",q="Q",r="R",
                      s="S",t="T",u="U",v="V",w="W",x="X",y="Y",z="Z")
      if (length(u.l) == 0) u.l<-sstr
      #print(c(sstr,u.l,upper.case))
      upper.case[is]<-paste(upper.case[is],u.l,sep="")
    }
  }
  if (lfac) {
    upper.case<-factor(upper.case)
  }
  upper.case
}
  
