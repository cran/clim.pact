# Estimates the modulo of two numbers: mod(y1,y2);
# Eg. mod(3,12)=3; mod(20,10)=0; mod(13,6)=1;
# R.E. Benestad, DNMI, 04.01.2001
#
mod <- function(x, y) {
  x1<-trunc( trunc(x/y)*y )
  z<-trunc(x)-x1
  z
  }
