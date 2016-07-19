distancia=function(A,B){
  "retorna una tabla con las distancias entre A y B fila vs columna en ese orden"
  x1=A$Long   
  y1=A$Lat
  x2=B$Long
  y2=B$Lat
  
  d1=as.matrix(cbind(x1,y1))
  d2=as.matrix(cbind(x2,y2))
  A=spDists(d1,d2,longlat = T)
  A
}

validacionTablas=function(rasterList,umbral,obs,auto,inter=1,direc){
  
  # Las estaciones automaticas a usar tienen que estar a una distancia umbral de las estaciones usadas
  # para la generacion del pisco por ello se compara con una matriz de distancia
  dias=0 #dias que se cumple la condicion de disntacia
  fechas=seq(as.Date('2014-01-01'),as.Date('2015-12-31'), "day")
  a=c(0,0,0,2,5,10,20,50)
  b=c(0,0,2,5,10,20,0,0)
  for (i in 4:733){
    piscopluvio=data.frame(Long=obs$x,Lat=obs$y,pp=obs[[i]])
    piscopluvio=na.omit(piscopluvio)
    
    automapluvio=data.frame(Long=auto$x,Lat=auto$y,pp=auto[[i]])
    automapluvio=na.omit(automapluvio)
    if (inter==1 || inter==7 || inter == 8){
      automapluvio=subset(automapluvio,pp>=a[inter])
    }else if (inter==2 ){
      automapluvio=subset(automapluvio,pp>a[inter])
    }else{
      automapluvio=subset(automapluvio,(pp>=a[inter] & pp<b[inter]))
    }
    N=dim(automapluvio)
    if (N[1]==0){
      print("No hay estaciones Automaticas")
    }else{
      
      A=distancia(automapluvio,piscopluvio)
      condicion=A>=umbral
      m=dim(A)
      selec=rep(T,m[1])#1:m[1]
      for (k in 1:m[1]){
        # si selec es True esa estacion no se usa en la validacion
        a=condicion[k,]
        selec[k]=(length(a[a==F])>0)
      }
      automapluvio$selec=selec
      automapluvio$distancia=apply(A,1,min)
      autousar=subset(automapluvio,selec==F) #subset que se usara para validar el producto
      autousar$selec=NULL
      puntos=autousar
      NN=dim(puntos)
      if (NN[1]==0){
        print('no hay estaciones')
      }else{
        coordinates(puntos)=~Long+Lat
        #valores extraidos del producto pisco para las posiciones donde las automaticas se ubican
        n=length(puntos)
        PPestimado=1:n
        if (n==0){
          print("No hay estaciones automaticas a la distancia correcta")
        }else{
          valores=extract(raster(rasterList[i-3]),puntos)
          for (j in 1:n){
            PPestimado[j]=valores[[j]]
          }
          resultado=cbind(autousar,PPestimado)
          resultado$error=abs(resultado$pp-resultado$PPestimado)
          resultado$landa=abs(resultado$pp-resultado$PPestimado)/resultado$pp
          write.table(resultado,paste0(paste0(dir,"/intervalo",as.character(inter)),'/',fechas[i-3],'.csv'),row.names = F,sep = ',')
          dias=dias+dim(resultado)[1]
        }
      }
    }
  }
  dias
}


cuadro=function(dir,inter){
  setwd(dir)
  errores=list.files(pattern = "\\.csv",recursive = T,include.dirs = T)
  
  
  
  rango1=0
  rango2=0
  rango3=0
  rango4=0
  rango5=0
  rango6=0
  
  
  
  for (j in 1:length(errores)){
    tablita=read.csv(errores[j])
    rango1 = rango1+ length(tablita$error[tablita$error==0])
    rango2 = rango2+length(tablita$error[tablita$error>0 & tablita$error<=0.5])
    rango3 = rango3+ length(tablita$error[tablita$error>0.5 & tablita$error<=2])
    rango4 = rango4+ length(tablita$error[tablita$error>2 & tablita$error<=5])
    rango5 = rango5+ length(tablita$error[tablita$error>5 & tablita$error<=10])
    rango6 = rango6+ length(tablita$error[tablita$error>10])
  }
  
  suma=rango1+rango2+rango3+rango4+rango5+rango6
  
  #numero de observaciones dias=66422
  
  resultados=data.frame(dias=suma,e1=rango1/suma,e2=rango2/suma,e3=rango3/suma,e4=rango4/suma,e5=rango5/suma,e6=rango6/suma)
  
  write.table(resultados,paste0('tabla',as.character(inter),'.csv'),sep = ',',row.names = F)
  
}