# R.E. Benestad, met.no, Oslo, Norway 04.06.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# NORDKLIMstations.

avail.elem <- function() {

ele <- c(101,111,112,113,121,122,123,401,601,602,701,801,911)
ele.c <- c("TAM","TAX","Th","Thd","TAN","Tl","Tld","SLP","RR","RRX",
           "DSC","CLOUD","SDM")
nam <- c('mean T(2m)','mean maximum T(2m)','highest maximum T(2m)',
         'day of Th date Thd','mean minimum T(2m)','lowest minimum T(2m)',
         'day of Tl date Tld','mean SLP','monthly accum. precip.',
         'maximum precip.',
         'Number of days with snow cover (> 50% covered) days dsc',
         'Mean cloud cover % N',
         'mean snow depth')


avail.elem <- list(data.set="Nordklim",ele=ele,name=nam)
avail.elem
}
