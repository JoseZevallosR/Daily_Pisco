
library(shiny)
library(leaflet)
library(raster)
library(tcltk)
path0='C:/Control'
source(paste0(path0,"/funcionesPisco.R"))
source(paste0(path0,"/download_GPM.R"))
source(paste0(path0,"/trmmfinal.R"))

source(paste0(path0,"/PiscoAutomatico.R"))

server <- function(input, output, clientData,session) {

  observeEvent(input$merge,{
    #path=tk_choose.dir()
    interpolacion(path0)
  })
  observeEvent(input$des, {
    path=tk_choose.dir()
    electionR <- input$dataset1
    
    c_num1 <- input$date
    c_num2 <- input$date2
    
    if (electionR=="RT3B42"){
      DOWNF(inicio = as.character(c_num1)  ,fin =as.character(c_num2),rtg=path)
    }else{
      downloadGPM(fileLocation = paste0(path,'/'),
                  user="aybar1994@gmail.com",
                  password="aybar1994@gmail.com",
                  dates=c(fil(c_num1),fil(c_num2)),
                  typeProduct = "early",
                  maxspeeddownload= 504800, # 504,800 Bytes/second (500KB/s)
                  numconnections='3',
                  band=6,
                  BBlonMin = -86,
                  BBlonMax = -66,
                  BBlatMin = -19.25,
                  BBlatMax = 1.25)        
    }
    
  })
  #})
  observeEvent(input$ppt,{
    output$mymap <- renderLeaflet({
      
      r=raster(choose.files())
      pal <- colorNumeric(c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                            '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'), values(r))
      leaflet() %>% addTiles() %>%
        
        
        addRasterImage(r, colors = pal, opacity = 0.4,project = F) %>%
        addLegend(pal = pal, values = values(r),
                  title = "lenyenda")
  })

    
  })
  
  observeEvent(input$extraccion,{
    #archivo de las previsiones
    
    archivo=read.table(choose.files(),header = F)
    #archivo del shape file
    #shape=choose.files()
    #abre el shape file
    puntos=shapefile(choose.files())
    #inicializa la grilla para el peru resolucion 0.25
    
    ###realiza la interpolacion
    setwd(paste0(path0,'/prevision'))
    ras=interpolacions(archivo)
    #guarda los puntos extraidos del raster
    
    resultado=extraccion(puntos,ras)
    names(resultado)=c('nombre','x','y','z','h24','h48','h72')
    
    write.table(resultado,paste0(Sys.Date(),'.csv'),sep=',',row.names = F)
  })
  
  
}
