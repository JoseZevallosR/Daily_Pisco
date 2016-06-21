################################################################
#Programa para descargar datos TRMM y acumularlos a paso Diario#
# Autor: Cesar Luis Aybar Camacho, SENAMHI
# Uso:   libre
################################################################
#'@description This function need to declare the start and ends of the period of interest, additionally 
#'you need to declare the working directory in "rtg".

DOWNF<-function(inicio="2015-02-28",fin="2015-03-08",product="TRMM",rtg="C:/Users/USUARIO//Desktop/1/"){ 
  if (!require(RCurl)) stop("Package RCurl is not installed")
  if (!require(raster)) stop("Package raster is not installed")
  if (!require(rgdal)) stop("Package rgdal is not installed")
  if (!require(gdalUtils)) stop("Package gdalUtils is not installed")
if(product=="TRMM"){ 
#---------------------------------------
fechadeinicio<-as.Date(inicio)
fechadefin<-as.Date(fin)
#---------------------------------------
rutadeguardado<-rtg  #define working directory

#-------------------------------------------------------
#Create a list with the dates to save
#-------------------------------------------------------
sq<-seq(from = fechadeinicio,to = fechadefin,by = "day")
levels<-levels(factor(substr(sq,1,7)))
year<-substr(levels,1,4)
month<-substr(levels,6,7)
url<-"ftp://disc2.nascom.nasa.gov/data/s4pa/TRMM_RT/TRMM_3B42RT.007"
howfile<-function(i){
myURL <- paste(url,"/",year[i],"/","",month[i],"/",sep="")
filenames <- getURL(myURL, ftp.use.epsv = FALSE, ftplistonly=TRUE, crlf=TRUE)
filePaths <- paste(myURL, strsplit(filenames, "\r*\n")[[1]], sep="")
#La siguiente linea selecciona solo  los productos .bin
selectedfilePaths <- filePaths[grep(filePaths,pattern=paste("\\.bin$"))] }
k<-mapply(howfile,i=1:length(year))  
baseF<-unlist(k)
baseF<-sort(baseF)
#-------------------------
#Download day for day
#-------------------------
DownTRMM<-function(i){
day1<-baseF[grep(baseF,pattern = paste0("\\.",substr(sq[i],1,4),substr(sq[i],6,7),substr(sq[i],9,10)))]
day2<-baseF[grep(baseF,pattern = paste0("\\.",substr(sq[i]+1,1,4),substr(sq[i]+1,6,7),substr(sq[i]+1,9,10)))]
dayF<-c(day1,day2)
selectedfilePaths<-dayF[5:12]
#Esta linea empieza la descarga 
dir.create(paste0(rutadeguardado,"/",sq[i]))
dir.create(paste0(rutadeguardado,"/",sq[i],"/binaries"))
dir.create(paste0(rutadeguardado,"/",sq[i],"/3h"))
setwd(paste0(rutadeguardado,"/",sq[i],"/binaries"))
mapply(download.file, selectedfilePaths, basename(selectedfilePaths))

#--------------------------------------------------------------------------------------------
# Script extracted from "Read_trmm_3B42RT_1.R" create by Carlos Antonio Fernandez-Palomino
# E-mail: cafpxl@gmail.com
#--------------------------------------------------------------------------------------------

fili = Sys.glob("*.bin")
fili_prefix = "TRMM.3B42RT"#
nlat   = 480
mlon   = 1440
lhead  = 2880     # characters (bytes)
prcScale = 0.01    #Scale factor to prec 
len<-length(fili)

for (j in 1:len){
  # extract data: typical file name: 3B42RT.2007062518.bin
  #                                  01234567890123456
  yyyy           = substr(fili[j], 8,11) 
  mm             = substr(fili[j],12,13) 
  dd             = substr(fili[j],14,15) 
  hh             = substr(fili[j],16,17) 
  yyyymmddhh     = substr(fili[j], 8,17) 
  print(paste(fili[j],"  nf=","for","  yyyy=",yyyy,"  mm=",mm,"  dd=",dd,"  hh=",hh))
  
  
  #***************************************************************
  # Read precip (short==integer)# skip header# read type short# reshape & scale
  # flag_value=-31999   (type short)
  #      integer precipitation and random error fields are clipped to 
  #      [-31998,31998] to prevent duplication of the missing value 
  #      (at the negative end)
  #***************************************************************
  
  lhead2 = lhead/2  # header length as type 'short' 
  nread_bin = lhead2+nlat*mlon# header length + prec length
  #read binary data
  works = readBin(paste(rutadeguardado,'/',sq[i],"/binaries/",fili[j],sep=""), n=nread_bin,what="integer",endian = "big",size=2,signed = TRUE)
  #works = readBin("3B42RT.2001010200.7R2.bin", n=691200,what="integer",endian = "big",size=2,signed = TRUE)
  
  #head(works,1450)#verify
  #hist(works)#histogram
  #min(works,na.rm = T)
  works[works==-31999] <- NA#missing values
  #hist(works)#histogram
  #only prec value
  prec<-works[(lhead2+1):nread_bin]*prcScale# Scale correction 
  #hist(prec);min(prec,na.rm = T);max(prec,na.rm = T)
  
  #***************************************************************
  # Eliminate negative precipitation values for graphics
  #***************************************************************
  prec[prec<0] <- NA#corregiendo valores negativos
  #hist(prec);min(prec,na.rm = T);max(prec,na.rm = T)
  
  #***************************************************************
  # create raster, specify extent and projection
  #***************************************************************
  
  baseraster <- raster(extent(-180,180,-60,60))
  res(baseraster) <- 0.25
  projection(baseraster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  #filling row acording:
  #first_box_center=59.875N,0.125E 
  #second_box_center=59.875N,0.375E 
  #last_box_center=59.875S,359.875E
  mtrx = matrix(prec, nrow=nlat, ncol=mlon, byrow=T)
  
  prec_ok=cbind(mtrx[,721:1440],mtrx[,1:720])
  baseraster[] <- prec_ok#longitude from -180 to 180, validated how in ncl
  #baseraster[] <- mtrx#de 0 a 360 validated how in ncl language
  
  #plot(baseraster)#displaying
  
  #study area
  area <- crop(baseraster,extent(-86,-66,-19.25,1.25))#PERU PERSIAN 0.04 res extent(-86,-66,-19.20,1.20), VILCANOTA extent(-73.4,-70.4,-15.2,-12.4) para vilcanota
  #plot(area)
  #***************************************************************
  # Write in format Geotiff
  #***************************************************************
  writeRaster(area, filename=paste0(rutadeguardado,"/",sq[i],"/3h","/",substr(fili,1,19),".tif")[j], overwrite=TRUE)
}

setwd(paste0(rutadeguardado,'/',sq[i],"/3h/"))
#Prec daily
fili2 = Sys.glob("*.tif")
dialy<-sum(stack(fili2))
setwd(paste0(rutadeguardado,'/',sq[i]))
writeRaster(dialy,paste("TRMM.3B42RT.",substr(fili2[i],8,15),".tif",sep=""))}
k<-mapply(DownTRMM,i=1:length(sq))  
}}
#EXAMPLE


