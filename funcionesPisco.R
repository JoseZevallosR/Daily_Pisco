library(zoo)
library(gstat)
library(automap)
library(sp)
library(rgdal)
library(hydroGOF)
library(raster)
require(gdalUtils)
require(RCurl)
library(XML)
library(hydroGOF)


#path="/home/senamhi-01/Documentos/data_products/Entregable"
######################################################################################################
descargas=function(mes,ano){ 
  library(RSelenium)
  #abriendo el server
  
  RSelenium::startServer()
  checkForServer()
  
  remDr <- remoteDriver(browserName = "chrome")
  remDr$open()
  # abre la pagina
  url<-"http://www.senamhi.gob.pe/site/lvera/red_rpm2.php"
  remDr$navigate(url)
  meses=c('Enero','febrero','marzo','abril','mayo','junio','julio','agosto','setiembre','octubre','noviembre','diciembre')
  ######
  j=paste("//*/option[@value = ",as.character(mes),']',sep = '')
  mesE <- remDr$findElement(using = 'xpath',j )
  mesE$clickElement()
  
  ####
  k=paste("//*/option[@value = ",as.character(ano),']',sep = '')
  anoE <-remDr$findElement(using = 'xpath', k)
  anoE$clickElement()
  
  # identify search button and click
  searchID<-'//*[@name="test2"]'
  webElem<-remDr$findElement(value = searchID)
  webElem$clickElement()
  
  # identify the table
  doc <- htmlParse(remDr$getPageSource()[[1]])
  tabla=readHTMLTable(doc)[[3]]
  nombre=paste(meses[mes],as.character(ano),'.csv',sep='')
  
  #acabado final
  fecha=seq(as.Date("2011/01/01"),as.Date("2016/12/31"),'day')
  write.table(tabla,nombre,sep=',',row.names = F,col.names = F)
  remDr$close()
}


tablas=function(){
  meses=c('Enero','febrero','marzo','abril','mayo','junio','julio','agosto','setiembre','octubre','noviembre','diciembre')
  file=paste(meses[mes],as.character(ano),'.csv',sep='')
  tabla<-read.csv(file)
  PP24<-tabla[tabla$VAR=="PP24",]
  dir=path
  setwd(dir)
  write.csv(PP24,paste(mes,ano,"_PP",".csv",sep=""))
}

vozdata=function(){
  ur='http://www.senamhi.gob.pe/site/lvera/reporte_diario.php'  
  tabla=readHTMLTable(ur,header = T)[[2]]
  data=tabla[2:251,3:5]
  V3=as.vector(data$V3)
  V4=as.vector(data$V4)
  data$V3=as.numeric(V4)
  data$V4=as.numeric(V3)
  data$V5=as.numeric(data$V5)
  row.names(data)=NULL
  names(data)=c('x','y','obs')
  data=na.omit(data)
  data$obs=log1p(data$obs)
  data
}

#Calcular el pbias entre un raster(satelite) y un spatialpointsdataframe (estaciones)
pbiasSAT<-function(gauge=dperu2,sat=rSAT){
  ext<-extract(sat,gauge,cellnumber=F,sp=T)
  names(ext)<-c("gauge","sat")
  #r<-cor(as.numeric(ext[[5]]),ext[[6]]*+10, use = "pairwise.complete.obs")
  pvalue<-cor.test(as.numeric(ext[[1]]),ext[[2]], method = "pearson")[3]
  return(pvalue$p.value)
}

#############################################
#########PROMEDIAR ESTACIONES DOBLES#########
#############################################

mean.doble.Station=function (gauge=dperu2, sat=rSAT, longlat = TRUE) 
{ 
  spdf<-gauge
  ts<-gauge@data
  point<-rasterToPoints(sat)
  point<-as.data.frame(point)
  colnames(point)<-c("x","y","sat")
  coordinates(point)<-~x+y
  projection(point)<-projection(sat)
  loc<-as.numeric()
  for (i in 1:length(gauge))     loc[i] <- which.min(spDists(point, gauge[i,], longlat))
  locsave<-loc  
  ts_colocated <- matrix(ncol = length(unique(loc)), nrow = 1)
  for (t in 1:nrow(t(ts))) ts_colocated[t, ] <- tapply(ts[, t], loc, mean, na.rm = TRUE)
  ts_colocated[!is.finite(ts_colocated)] <- NA
  length_isna <- numeric()
  for (i in 1:ncol(ts_colocated)) length_isna[i] <- length(which(is.na(ts_colocated[, i])))
  remove <- which(length_isna == nrow(ts_colocated))
  if (length(remove > 1)) 
    ts_colocated <- ts_colocated[, -remove]
  
  ts_colocated <- as.zoo(ts_colocated)
  index(ts_colocated) <- index(t(ts))
  names(ts_colocated) <- 1:ncol(ts_colocated)
  gauge <- point[unique(loc)[order(unique(loc))], ]
  if (length(remove > 1)) 
    gauge <- gauge[-remove, ]
  gauge$gauge <- as.numeric(ts_colocated)
  gauge2<-gauge[2]
  names(gauge2)<-names(spdf)
  gauge<-gauge2
  return(gauge)
}

###############################################
################CALCULAR VARIOGRAMA############
###############################################

variogrm<-function(gauge=dperu2,sat=rSAT){ 
  ext<-extract(sat,gauge,cellnumber=F,sp=T)
  names(ext)<-c("gauge","sat")
  vm.fit<- afvmod(gauge~sat,ext,fix.values = c(0,NA,NA),cressie = TRUE)
  return(vm.fit)
}

KED<-function(gauge=dperu2,sat=rSAT){ 
  ext<-extract(sat,gauge,cellnumber=F,sp=T)
  names(ext)<-c("gauge","sat")
  vm.fit<- afvmod(gauge~sat,ext,fix.values = c(0,NA,NA),cressie = TRUE)
  #save(vm.fit,file=paste('/home/senamhi-cesar/Escritorio/newDailyInterpolation/KED/variogrm/Variograma',substr(names(sat),13,16),".",substr(names(sat),18,19),".",substr(names(sat),21,23),".rda",sep=""))
  #png(paste('/home/senamhi-cesar/Escritorio/newDailyInterpolation/KED/plotFitted//Variograma',substr(names(sat),13,16),".",substr(names(sat),18,19),".",substr(names(sat),21,23),".png",sep=""))
  message("Grabando variograma")
  print(plot(vm.fit))
  #Define grid
  point<-rasterToPoints(sat)
  point<-as.data.frame(point)
  colnames(point)<-c("x","y","sat")
  coordinates(point)<-~x+y
  projection(point)<-projection(sat)
  Zs<-krige(gauge~sat, locations = ext, newdata = point, model = vm.fit$var_model)
  map   		<- as(Zs[1],"SpatialPixelsDataFrame")
  gridded(map) <- TRUE
  mapa <- raster(map)
  mapa<-expm1(mapa)
  mapa[mapa<=0]=0
  writeRaster(mapa,filename  = paste0('PISCO.PREC.V1.',substr(names(sat),13,16),".",substr(names(sat),17,18),".",substr(names(sat),19,20),".tif"),overwrite=TRUE)
}


#####################################
################IDW##################
#####################################


IDW<-function(gauge=dperu2,sat=rSAT){ 
  #dperu<-cbind(dperu[1:6],dperu[7678:length(colnames(dperu))])
  ext<-extract(sat,gauge,cellnumber=F,sp=T)
  point<-rasterToPoints(sat)
  point<-as.data.frame(point)
  colnames(point)<-c("x","y","sat")
  coordinates(point)<-~x+y
  projection(point)<-projection(sat)
  names(ext)<-c("gauge","sat")
  Zs<-idw(expm1(gauge)~1, locations = ext, newdata = point)
  map   		<- as(Zs[1],"SpatialPixelsDataFrame")
  gridded(map) <- TRUE
  mapa <- raster(map)
  mapa[mapa<=0]=0
  writeRaster(mapa,filename  = paste0('PISCO.PREC.V1.',substr(names(sat),13,16),".",substr(names(sat),17,18),".",substr(names(sat),19,20),".tif"),overwrite=TRUE)
  
}

#####################################
########Regression-IDW###############
#####################################



RegressionIDW<-function(gauge=dperu2,sat=rSAT){ 
  ext<-extract(sat,gauge,cellnumber=F,sp=T)
  names(ext)<-c("gauge","sat")
  station<-gauge
  llm<-lm(ext[[1]]~ext[[2]])
  data<-match(seq(1:946),as.numeric(names(llm$residuals)))
  joi<-as.numeric()
  station$residuals<-llm$residuals
  grd.pts <- SpatialPixels(SpatialPoints((sat)))
  grd <- as(grd.pts, "SpatialGrid")
  projection(grd)<-projection(station)
  idwError<-idw(residuals~1,station,grd,idp=2)
  idwError<-idwError["var1.pred"]
  
  #guardando el tiff
  gridded(idwError) <- TRUE
  mapa <- raster(idwError)
  OBSp<-sat*llm$coefficients[2]+llm$coefficients[1]
  Ridw<-OBSp+mapa
  Ridw[Ridw<0]=0
  Ridw=expm1(Ridw)
  writeRaster(Ridw,filename  = paste0('PISCO.PREC.V1.',substr(names(sat),13,16),".",substr(names(sat),17,18),".",substr(names(sat),19,20),".tif"),overwrite=TRUE)
}


###################################################
################AUTOFIT VARIOGRAM##################
###################################################

afvmod = function(formula, input_data, model = c("Sph", "Exp", "Gau", "Ste"),
                  kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA),
                  verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA), 
                  miscFitOptions = list(), fit.method = 7, ...)
  # This function automatically fits a variogram to input_data
{
  # Check for anisotropy parameters
  if('alpha' %in% names(list(...))) warning('Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.')
  
  # Take the misc fit options and overwrite the defaults by the user specified ones
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5, 
                                orig.behavior = TRUE, 
                                num.bins = NA, equal.width.bins = FALSE,
                                equal.np.bins = FALSE)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  
  # Create boundaries
  longlat = !is.projected(input_data)
  if(is.na(longlat)) longlat = FALSE
  diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1,2] * 0.35 # 0.35 times the length of the central axis through the area
  
  ### BEGIN MODIFICATIONS ###
  
  if(miscFitOptions[["orig.behavior"]]){
    if(verbose) cat ("Boundaries as defined by original autofitVariogram...\n\n")
    # compute boundaries the old way
    boundaries = c(5,10,15,20,25,30,50,70,90,100,130,160) # Boundaries for the bins in km
    # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
    # request by Jon Skoien
    if(miscFitOptions[["merge.small.bins"]]) {
      # bin the old way
      if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while(TRUE) {
        if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
        boundaries = boundaries[2:length(boundaries)]                 
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
        } else {
          experimental_variogram = variogram(g, boundaries = boundaries, ...)
        }
      } 
    }
    ### equal-width bins (approximately, variogram does its own binning too) ###
  } else if(miscFitOptions[["equal.width.bins"]]){
    if(verbose) cat("Using equal-width bins...\n")
    if('width' %in% names(list(...))) stop('Cannot pass width when equal.width.bins = TRUE. Supply "init.width" in miscFitOptions instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # user must supply either width or num.bins
    if(!'init.width' %in% names(miscFitOptions)){
      if(is.na(miscFitOptions[['num.bins']])) stop('when equal.width.bins = TRUE, user must supply either init.width or num.bins as well.')
      width <- diagonal/miscFitOptions[['num.bins']]
      if(verbose) cat("initial width not provided. Calculating using num.bins.\n")
    } else {
      width <- miscFitOptions[['init.width']]
      if(verbose) cat("initial width provided.\n")
    }
    # get the experimental variogram
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, width = width, ...)
    } else {
      experimental_variogram = variogram(g, width = width, ...)
    }
    # merge small bins if requested
    if(miscFitOptions[['merge.small.bins']]){
      if(verbose) cat("Checking if any bins have less than ", miscFitOptions[["min.np.bin"]], " points, merging bins when necessary...\n")
      iter <- 0
      maxiter <- 1000
      while(TRUE){
        if(!any(experimental_variogram$np < miscFitOptions[["min.np.bin"]])) break  		
        # increase width by 10% and try again
        width <- width*1.1
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, width = width, ...)
        } else {
          experimental_variogram = variogram(g, width = width, ...)
        }
        iter <- iter + 1
        if(iter > maxiter){
          cat('maximum number of interations reached. Try decreasing min.np.bin or init.width.\n\n')
          break
        }
      }
    }
    ### equal observation count bins ###
  } else if(miscFitOptions[["equal.np.bins"]]){
    if(verbose) cat("Using bins of equal observation counts...\n")
    if('boundaries' %in% names(list(...))) stop('Cannot pass boundaries when equal.np.bins is TRUE. Pass num.bins or min.np.bin instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # get a sorted list of distances 
    dists <- sort(spDists(input_data))
    # apply the cutoff
    dists <- dists[dists < diagonal & dists > 0]
    # split the data into bins based on number of observations
    if(is.na(miscFitOptions[['num.bins']])){
      # compute number of bins based on the minimum number of observations per bin
      miscFitOptions[['num.bins']] <- floor(0.5*length(dists)/miscFitOptions[['min.np.bin']])
      if(verbose) cat("num.bins not supplied. Setting num.bins =", miscFitOptions[['num.bins']], 'based on min.np.bin.\n')
    }
    cat("checking bins, decreasing num.bins if necessary... \n")
    while(TRUE){
      # compute interval based on the number of bins
      interval <- length(dists)/miscFitOptions[['num.bins']]
      # define boundaries
      boundaries <- rep(NA, miscFitOptions[['num.bins']])
      for(i in 1:miscFitOptions[['num.bins']]){
        boundaries[i] <- dists[round(i*interval)]
      }
      if(length(boundaries == length(unique(boundaries)))) break
      # reduce number of bins
      miscFitOptions[['num.bins']] <- miscFitOptions[['num.bins']] - 1 
    }
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
  } else {
    # default behavior of variogram
    cat("No binning action specified in miscFitOptions.\n\n")
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, ...)
    }
  }
  ### END MODIFICATIONS ###
  
  # set initial values
  if(is.na(start_vals[1])) {  # Nugget
    initial_nugget = min(experimental_variogram$gamma)
  } else {
    initial_nugget = start_vals[1]
  }
  if(is.na(start_vals[2])) { # Range
    initial_range = 0.1 * diagonal   # 0.10 times the length of the central axis through the area
  } else {
    initial_range = start_vals[2]
  }
  if(is.na(start_vals[3])) { # Sill
    initial_sill = mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma)))
  } else {
    initial_sill = start_vals[3]
  }
  
  # Determine what should be automatically fitted and what should be fixed
  # Nugget
  if(!is.na(fix.values[1]))
  {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  } else
    fit_nugget = TRUE
  
  # Range
  if(!is.na(fix.values[2]))
  {
    fit_range = FALSE
    initial_range = fix.values[2]
  } else
    fit_range = TRUE
  
  # Partial sill
  if(!is.na(fix.values[3]))
  {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  } else
    fit_sill = TRUE
  
  getModel = function(psill, model, range, kappa, nugget, fit_range, fit_sill, fit_nugget, fit.method, verbose)
  {
    if(verbose) debug.level = 1 else debug.level = 0
    if(model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if(is.na(start_vals[1])) nugget = 0
      if(is.na(start_vals[2])) range = 1    # If a power mode, range == 1 is a better start value
      if(is.na(start_vals[3])) sill = 1
    }
    if(fit.method==5){
      #fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
      #	obj = try(fit.variogram.reml(experimental_variogram,
      #					        model=vgm(psill=psill, model=model, range=range,
      #							nugget=nugget,kappa = kappa),
      #			  fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
      #			  debug.level = 0, fit.method=fit.method), TRUE)		
    }
    else if(fit.method==8){
      #fit.variogram.gls(formula, data, model, maxiter = 30,
      #eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
      #plot = FALSE
      #obj = try()
    }
    else{
      obj = try(fit.variogram(experimental_variogram,
                              model=vgm(psill=psill, model=model, range=range,
                                        nugget=nugget,kappa = kappa),
                              fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
                              debug.level = 0, fit.method=fit.method), TRUE)
    }
    if("try-error" %in% class(obj)) {
      #print(traceback())
      warning("An error has occured during variogram fitting. Used:\n", 
              "\tnugget:\t", nugget, 
              "\n\tmodel:\t", model, 
              "\n\tpsill:\t", psill,
              "\n\trange:\t", range,
              "\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
      return(NULL)
    } else return(obj)
  }
  
  
  # Automatically testing different models, the one with the smallest sums-of-squares is chosen
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  
  for(m in test_models) {
    if(m != "Mat" && m != "Ste") {        # If not Matern and not Stein
      model_fit = getModel(initial_sill - initial_nugget, m, initial_range, kappa = 0, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, fit.method=fit.method)
      if(!is.null(model_fit)) {       # skip models that failed
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
      counter = counter + 1
    } else {                 # Else loop also over kappa values
      for(k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget, m, initial_range, k, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, fit.method=fit.method)
        if(!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
        counter = counter + 1
      }
    }
  }
  
  # Check for negative values in sill or range coming from fit.variogram
  # and NULL values in vgm_list, and remove those with a warning
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill, v$range) < 0) | is.null(v))
  if(any(strange_entries)) {
    if(verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  
  if(verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame("Tested models" = sapply(vgm_list, function(x) as.character(x[2,1])), 
                        kappa = sapply(vgm_list, function(x) as.character(x[2,4])), 
                        "SSerror" = SSerr_list)
    tested = tested[order(tested$SSerror),]
    print(tested)
  }
  
  result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]], sserr = min(SSerr_list))
  class(result) = c("autofitVariogram","list")    
  
  return(result)
}

fil=function(a){
  a=paste0(substr(a,1,4),'.',substr(a,6,7),'.',substr(a,9,10))
  a
}

tablaExtraccion=function(ras,cuenca){
  datos=as.data.frame(rasterToPoints(ras))
  datos$X=datos$x
  datos$Y=datos$y
  coordinates(datos)=~x+y
  proj4string(datos)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0  '
  Estaciones=over(cuenca,datos,returnList = T)
  CUENCA=Estaciones$`0`
  tabla=data.frame(x=CUENCA$X,y=CUENCA$Y,pp=CUENCA[[1]])
  write.table(tabla,'datos.csv',sep = ',',row.names = F)
}

grilla=function(resX,resY){
  pixelsx=(86-66)/resX
  pixelsy=(19.25-1.25)/resY
  geog.grd <- expand.grid(x=seq(floor(-86),
                                ceiling(-66),
                                length.out=pixelsx),
                          y=seq(floor(-19.25),
                                ceiling(1.25),
                                length.out=pixelsy))
  grd.pts <- SpatialPixels(SpatialPoints((geog.grd)))
  grd <- as(grd.pts, "SpatialGrid")
  proj4string(grd)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0  '
  grd
}

interpolacions=function(datapoints){
  grd=grilla(0.15,0.15)
  names(datapoints)=c('y','x','h24','h48','h72')
  coordinates(datapoints)=~x+y
  proj4string(datapoints)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0  '
  
  fecha=Sys.Date()
  
  fechas=seq(fecha, fecha+2, "day")
  result=list(0) #lista que alamacenara los rasters
  
  for(i in 1:3 ){
    prueba<- afvmod(datapoints[[i]]~1,datapoints,fix.values = c(0,NA,NA),cressie = TRUE)
    #comprueba que el variograma este bueno
    if (attr(prueba$var_model,"singular")==FALSE & prueba$var_model$range[2]<150 & prueba$var_model$range[2]>7&prueba$sserr>0.0001){
      Zs<-krige(datapoints[[i]]~1, locations = datapoints, newdata = grd, model = prueba$var_model)
    }else{
      Zs <- idw(datapoints[[i]] ~ 1, datapoints, grd, idp=2)
    }  
    result[i]=raster(Zs)
    writeGDAL(Zs,paste(fechas[i],'.tif',sep = ''))
  }
  result
}

extraccion=function(puntos,ras){
  #shape : puntos de estaciones
  #ras  : lista de rasters
  #extraer los valores de raster a los puntos
  ext=lapply(ras,extract,puntos)
  m=length(ras)
  n=length(puntos)
  datos=matrix(data=0,nrow = n,ncol = m)
  for (i in 1:m){
    datos[,i]=ext[[i]]
  }
  #el shape tiene que tener en su tabla de atributos nombre,x,y,z
  resultado=data.frame(nombre=puntos[[1]],x=puntos[[2]],y=puntos[[3]],z=puntos[[4]])
  resultado=cbind(resultado,datos)
  resultado
}


tablas=function(obs1,prod){
  estadisticos=matrix(data=0,nrow = dim(obs1)[1],ncol = 4)
  for (i in 1:dim(obs1)[1]){
    nn=dim(obs1)[2]
    sim=as.numeric(prod[i,5:nn])
    obs=as.numeric(obs1[i,5:nn])
    f=gof(sim,obs,na.rm = T)
    estadisticos[i,]=c(f["RMSE",],f["r",],f["NSE",],f["PBIAS %",])
  }
  resultados=as.data.frame(estadisticos)
  resultados=cbind(obs1$NOMBRE,resultados)
  names(resultados)=c('Nombres','RMSE','r','NSE',"PBIAS")
  resultados
}