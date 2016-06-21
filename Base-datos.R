
path0='C:/Control'
source(paste0(path0,"/funcionesMinerve.R"))

setwd(paste0(path0,"/Base-Datos-Minerve"))

baseHistorica=read.csv('obrajillo.csv')
baseSimualacion=read.csv('obrajillo_prevision.csv')


#resolucon de la grilla para la prevision
grd=grilla(0.25,0.25)
#puntos de las estaciones manejados en la base de datos
setwd(paste0(path0,"/carlos-shape"))
puntos=shapefile('ESTACIONES_MODELO.shp')
#interpolacion de la prevision y extraccion de los valores para las estaciones de la base de datos
setwd(paste0(path0,"/prevision"))
prevision=extraccion(puntos,interpolacion(read.table('prevision.txt'),grd))

#creando la base actual
path_acumulado=paste0(path0,"/DIARIO")
path_gpm=paste0(path0,"/GPM")
path_baseDdatos=paste0(path0,"/Base-Datos-Minerve")

baseSim=baseActual(baseHistorica,prevision,puntos,path_acumulado,path_gpm,path_baseDdatos)
#DOWNFgpm(inicio="2016-02-28",fin="2016-02-28",fileLocation=path_gpm,fileAcum=path_acumulado)
