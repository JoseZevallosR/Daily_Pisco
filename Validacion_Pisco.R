library(raster)
distancia=function(A,B){
  "retorna una tabla con las distancias entre A y B fila vs columna en ese orden"
  x1=A$Long   
  y1=B$Lat
  x2=B$Long
  y2=B$Lat
  A=matrix(data=0,nrow =length(x1) ,ncol =length(x2) )
  for(i in 1:length(x1)){
    for(j in 1:length(x2)){
      A[i,j]=sqrt((x1[i]-x2[j])^2+(y1[i]-y2[j])^2)
    }
  }
  A
}

validacionTablas=function(rasterList,umbral,obs,auto){
  
  # Las estaciones automaticas a usar tienen que estar a una distancia umbral de las estaciones usadas
  # para la generacion del pisco por ello se compara con una matriz de distancia
  dias=0 #dias que se cumple la condicion de disntacia
  fechas=seq(as.Date('2014-01-01'),as.Date('2015-12-31'), "day")
  for (i in 4:733){
    piscopluvio=data.frame(Long=obs$x,Lat=obs$y,pp=obs[[i]])
    piscopluvio=na.omit(piscopluvio)
    
    automapluvio=data.frame(Long=auto$x,Lat=auto$y,pp=auto[[i]])
    automapluvio=na.omit(automapluvio)
    N=dim(automapluvio)
    if (N[1]==0){
      print("No hay estaciones Automaticas")
    }else{
      automapluvio=subset(automapluvio,(pp>=2 & pp<5))#indica el intervalo
      A=distancia(automapluvio,piscopluvio)
      condicion=A>=umbral
      m=dim(A)
      selec=rep(T,m[1])#1:m[1]
      for (k in 1:m[1]){
        # si selec es True esa estacion no se usa en la validacion
        a=condicion[k,]
        selec[k]=length(a[a==F])>1
      }
      automapluvio$selec=selec
      autousar=subset(automapluvio,selec==F) #subset que se usara para validar el producto
      autousar$selec=NULL
      puntos=autousar
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
        write.table(resultado,paste0('/home/jose/Documents/validacionPisco/tablas-pisco/intervalo4','/',fechas[i-3],'.csv'),row.names = F,sep = ',')
        dias=dias+dim(resultado)[1]
      }
    }
  }
  dias
}


validacionTablas=function(rasterList,umbral,obs,auto){
  
  # Las estaciones automaticas a usar tienen que estar a una distancia umbral de las estaciones usadas
  # para la generacion del pisco por ello se compara con una matriz de distancia
  dias=0 #dias que se cumple la condicion de disntacia
  fechas=seq(as.Date('2014-01-01'),as.Date('2015-12-31'), "day")
  for (i in 4:733){
    piscopluvio=data.frame(Long=obs$x,Lat=obs$y,pp=obs[[i]])
    piscopluvio=na.omit(piscopluvio)
    
    automapluvio=data.frame(Long=auto$x,Lat=auto$y,pp=auto[[i]])
    automapluvio=na.omit(automapluvio)
    automapluvio=subset(automapluvio,(pp>=0 & pp<2))#indica el intervalo
    #automapluvio=subset(automapluvio,pp>0)
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
        selec[k]=length(a[a==F])>1
      }
      automapluvio$selec=selec
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
          write.table(resultado,paste0('/home/jose/Documents/validacionPisco/tablas-pisco/intervalo3','/',fechas[i-3],'.csv'),row.names = F,sep = ',')
          dias=dias+dim(resultado)[1]
        }
      }
    }
  }
  dias
}