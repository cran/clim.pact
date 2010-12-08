narp.names <- c("Torshavn","Strond Kr.st.","Upernavik","Ilulissat Airport","Nuuk","Narsarsuaq",
                "Danmarkshavn","Ittoqqortoormiit","Tasiilaq","Stykkisholmur","Reykjavik",
                "Vestmannaeyar","Haell","Akureyri","Raufarhoefn","Teigarhorn","Ship M","Tromsoe",
                "Vardoe","Bjoernoeya","Hopen","Svalbard Airport","Svalbard Airport reconstructed",
                "Jan Mayen")
narp.lons <- c( -6.7667, -6.5833,-56.1667,-51.0667,-51.7500,-48.1667,-18.7667,-22.0000,-37.6333,
               -22.7333,-21.9000,-20.2833,-20.4167,-18.0833,-15.9500,-15.2000,  2.0000, 18.9333,
                31.0833, 19.0167, 25.0667, 15.4667, 15.4667, -8.6667)
narp.lats <- c( 62.0167, 62.2667, 72.7833, 69.2330, 64.1667, 61.2000, 76.7667, 70.4833, 65.6000,
                65.0833, 64.1333, 63.4000, 64.0667, 65.6833, 66.4500, 64.3000, 66.0000, 69.6500,
                70.3667, 74.5167, 76.5000, 78.2500, 78.2500, 70.9333)
narp.countries <- c(rep("Faeroes",2),rep("Greenland",7),rep("Iceland",7),rep("Norway",8))
narp.who.stnr <- c(6011,NA,4210,4221,4250,4270,4320,4339,4360,4013,4030,4048,4050,4063,4077,4092,NA,
              1026,1098,1028,1062,1008,NA,1001)
narp.stnr <- c(06011,33054,04210,04221,04250,04270,04320,04339,04360,04013,04030,04048,04050,
               04063,04077,04092,76900,90450,98550,99710,99720,99840,99841,99950)
narp.meta <- list(stnr=narp.stnr,names=narp.names,lons=narp.lons,lats=narp.lats,
                  countries=narp.countries,WMO.number=narp.who.stnr)
rm(narp.names,narp.lons,narp.lats,narp.countries,narp.who.stnr,narp.stnr)
gc(reset=TRUE)
