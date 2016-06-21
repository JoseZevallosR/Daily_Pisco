library(plyr)
library(snow)
library(raster)
library(rgdal)
library(doMC)
library(foreach)

setwd("//home/jose/Escritorio/TRMMdiario")
satlista=list.files(pattern = "\\.tif")
lim=1000#length(satlista)
x=stack(satlista[1:lim])


sumHarmonic=function(n){
  s=0
  scum=vector(mode="numeric", length=(n))
  for (i in n:1){
    s=s+i^-1
    scum[n-i+1]=s
  }
  scum
}



EP= function(x,n,weight){

  Harm=weight
  ep=0
  dim=length(x)
  for (N in n:1){
    if (is.na(x[dim-N+1])){
      ep=ep
    }else{
      ep=ep+x[dim-N+1]*Harm[n-N+1]
    }
  }
  ep
}
n=365
weight=sumHarmonic(n)## variable global

time.serieEP=function(x)
{
  n=365
  serie.EP=vector(mode="numeric", length=(length(x)-n))
  for (i in 1:(length(x)-n)){
    serie.EP[i]=EP(x[(n+i):(1+i)],n,weight)    
    
  }
  serie.EP
}

#21:33
EPTIF=calc(x,fun = time.serieEP)

plot(EPTIF)

MEP=calc(EPTIF,mean)

DEP=EPTIF-MEP
STDEP=calc(DEP,sd)
EDI=DEP/STDEP

plot(EDI)


setwd("/home/jose/Escritorio/EDI.Software_Linux_PGIFortran95/EDI_Mapas")
fecha=as.Date('2001-03-05')
for (i in 1:nlayers(EDI)){
  writeRaster(EDI[[i]],filename  =paste('EDI',fecha+i,'tif',sep = '.') ,overwrite=TRUE)
}







