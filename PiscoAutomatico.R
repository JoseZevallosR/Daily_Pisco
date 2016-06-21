library(raster)
library(zoo)

library(rgeos)

path0='C:/Control'
source(paste0(path0,"/funcionesPisco.R"))
source(paste0(path0,"/download_GPM.R"))
source(paste0(path0,"/trmmfinal.R"))
source(paste0(path0,"/accumIMERG.R"))



interpolacion=function(path){
  Fecha=Sys.Date()-1
  #DOWNF(inicio = as.character(Fecha),fin = as.character(Fecha),rtg=path)
  #para acumular descarga dos dias anteriores
  #downloadGPM(fileLocation = paste0(path,'/'),
  #            user="aybar1994@gmail.com",
  #            password="aybar1994@gmail.com",
  #            dates=c(fil(Fecha-1),fil(Fecha)),
  #            typeProduct = "early",
  #            maxspeeddownload= 504800, # 504,800 Bytes/second (500KB/s)
  #            numconnections='3',
  #            band=6,
  #            BBlonMin = -86,
  #            BBlonMax = -66,
  #            BBlatMin = -19.25,
  #            BBlatMax = 1.25)  
  #####IMERG Acumulado
  rSAT=AccumGPM( 
    InputLocation= paste0(path,'/GPM/') ,
    UTC="1200",
    TimeStep=  "Daily",
    paralelled=T,
    SO='Winbug',
    cores=7,
    OutputLocation=paste0(path,'/DIARIO/')
  )
  #data satelital
  
  
  setwd(paste0(path,'/DIARIO/'))
  lisSAT=list.files(pattern = "\\.tif")
  rSAT<-raster(lisSAT[length(lisSAT)])
  nameR<-names(rSAT)
  rSAT<-log1p(rSAT)
  names(rSAT)<-nameR
  #pluvio
  data=vozdata()
  coordinates(data)<-~x+y
  bubble(data,"obs")
  projection(data)<-projection(rSAT)
  pvalue<-pbiasSAT(gauge = data,sat = rSAT)
  dperu<- mean.doble.Station(gauge = data,sat = rSAT)
  #setwd(paste0(path,'/','Pisco'))
  setwd(paste0(path,'/'))
  if (file.exists('Pisco')){
    setwd(file.path(path, 'Pisco'))
  } else {
    dir.create(file.path(path, 'Pisco'))
    setwd(file.path(path, 'Pisco'))
    
  }
  if(pvalue>0.05||is.na(pvalue)) IDWsimple<-IDW(dperu,rSAT)
  if(pvalue<0.05) { 
    prueba<-try(variogrm(dperu,rSAT),silent = T)
    #png(paste0(path,'/KedTotal','/variogram.',substr(nameR,13,16),".",substr(nameR,17,18),".",substr(nameR,19,20),".png"))
    #print(plot(prueba))
    #dev.off()
    if (attr(prueba$var_model,"singular")==FALSE & prueba$var_model$range[2]<150 & prueba$var_model$range[2]>7&prueba$sserr>0.0001)  KED<-KED(dperu,rSAT)
    else  RIDW<-RegressionIDW(dperu,rSAT)
  }
}

