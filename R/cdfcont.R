# Obtains the name of variables in a netcdf File
#
# R.E. Benestad, 23.09.2003

cdfcont <- function(filename,path="") {

  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')

  system(paste(path,"ncdump -h  ",filename," > cdfcont.txt",sep=""),intern=T)
  cdfhead <- readLines("cdfcont.txt")
  cdfvars <- cdfhead[c(grep("float",lower.case(cdfhead)),
                       grep("short",lower.case(cdfhead)),
                       grep("double",lower.case(cdfhead)))]
  cdfdims <- cdfvars
  for (i in 1:length(cdfvars)) {
    i1 <- instring(" ",cdfvars[i])
    i2 <- instring("(",cdfvars[i])
    i3 <- instring(")",cdfvars[i])
    cdfdims[i] <- substr(cdfvars[i],i2[1]+1,i3[1]-1)
    cdfvars[i] <- substr(cdfvars[i],i1[1]+1,i2[1]-1)
  }
#  print(cdfvars)
#  print(cdfdims)
  torg <- cdfhead[grep("time_origin",lower.case(cdfhead))]
  if (length(torg)==0) {
    torg <- cdfhead[grep("since",lower.case(cdfhead))]
    t.org.pos <- regexpr("since",lower.case(torg))
    s<- instring('\"',torg)
    torg  <- substr(torg,t.org.pos+6,s[2]-1)
    dash <- instring("-",torg)
    spc <- instring(" ",torg)
    if (spc==0) spc <- nchar(torg)+1
    yy0 <- as.numeric(substr(torg,1,dash[1]-1))
    while (nchar(yy0) < 4) yy0 <- paste("0",yy0,sep="")
    mm0 <- as.numeric(substr(torg,dash[1]+1,dash[2]-1))
    dd0 <- as.numeric(substr(torg,dash[2]+1,spc[1]-1))
    while (nchar(dd0) < 2) dd0 <- paste("0",dd0,sep="")
    torg <- paste(dd0,cmon[mm0],yy0)
 #   print(paste("time.origin=",torg))
    if (is.na(dd0)) dd0  <- 15
  } else {
    s<- instring('\"',torg)
    if (length(s)==2) torg<-substr(torg,s[1]+1,s[2]-1)
  }
  tunit<- cdfhead[grep("time:unit",lower.case(cdfhead))]
  if (length(tunit)>0) {
     s<- instring('\"',tunit)
     if (length(s)==2) tunit<- strip(substr(tunit,s[1]+1,s[2]-1))
  } else tunit<-NULL

  offs <- cdfhead[grep("add_offset",lower.case(cdfhead))]
  if (length(offs)>0) {
    e <- instring('=',offs);  f <- regexpr("f;",offs) # f <- instring('f',offs)
    if (f <= 0) f <- regexpr("f ;",offs)
    if (f[1]>0) f<-f[length(f)] else f<-nchar(offs)
    yes <- (nchar(offs)>0) & (length(e)>0) & (length(f)>0)
    if (yes) offs<-substr(offs,e+1,f-1) else offs<-"0"
  } else offs <- "0"

  scal <- cdfhead[grep("scale_factor",lower.case(cdfhead))]
  if (length(offs)>0) {
    e <- instring('=',scal);  f <- regexpr("f;",scal) # f <- instring('f',scal)
    if (f <= 0) f <- regexpr("f ;",scal)
    if (f[1]>0) f<-f[length(f)] else f<-nchar(scal)
    yes <- (nchar(scal)>0) & (length(e)>0) & (length(f)>0)
    if (yes) scal<-substr(scal,e+1,f-1) else scal<-1
  } else scal<-1

  miss <- cdfhead[grep("missing_value",lower.case(cdfhead))]
  if (length(miss)>0) {
    if (length(miss) > 1) miss <- miss[1]
    e <- instring('=',miss);   f <- regexpr("f;",miss) #f <- instring('f',miss)
    if (f <= 0) f <- regexpr("f ;",miss)
    if (f[1]>0) f<-f[length(f)] else f<-nchar(miss)
    yes <- (nchar(miss)>0) & (length(e)>0) & (length(f)>0)
    if (yes) miss<-substr(miss,e+1,f-1) else miss<-NULL
  } else miss <- NA
  system("rm -f cdfcont.txt",intern=T)
  content <- list(vars=cdfvars,dims=cdfdims,time.origin=torg,time.unit=tunit,
                  add.offset=as.numeric(offs),scale.factor=as.numeric(scal),
       missing.value=as.numeric(miss))
  invisible(content)
}

