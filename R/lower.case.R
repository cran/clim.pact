# Returns the lower case of a string;
# Test:
# > lower.case("ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890") gives
# [1] "abcdefghijklmnopqrstuvwxyz1234567890"
# > lower.case("abcdefghijklemnoprstuvwxyz1234567890")
# [1] "abcdefghijklemnoprstuvwxyz1234567890"
# R.E. Benestad (REB), DNMI (met.no), 08.01.2002
#                REB:  03.05.2002 - modified to hande string arrays

lower.case <- function(u.case) {

  lfac<-FALSE                           # Set flag if we are dealing with a factor
                                        # object. Then the output is converted to
                                        # factor.
  
  if (is.factor(u.case)) { lfac <- TRUE }
  if (is.null(u.case)) return()
  str<-as.character(u.case)
#  print(str)
  lower.case<-str

  for (is in 1:length(str)) {
    nc<-nchar(str[is])
    lower.case[is]<-""
    for (ic in 1:nc) {
      sstr<-substr(str[is],ic,ic)
      if (sstr=="E") sstr<-"e"     # Fudge - E didn't work in switch..
      u.l<-switch(as.character(sstr),
                      A="a",B="b",C="c",D="d",F="f",G="g",H="h",I="i",
                      J="j",K="k",L="l",M="m",N="n",O="o",P="p",Q="q",R="r",
                      S="s",T="t",U="u",V="v",W="w",X="x",Y="y",Z="z")
      if (length(u.l) == 0) u.l<-sstr
#    print(c(sstr,u.l,lower.case))
      lower.case[is]<-paste(lower.case[is],u.l,sep="")
    }
  }
  if (lfac) {
    lower.case<-factor(lower.case)
  }
  lower.case
}
  
