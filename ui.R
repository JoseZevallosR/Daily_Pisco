library(shiny)
library(leaflet)
library(raster)

ui <- fluidPage(
  tabPanel("Interactive map",
           div(class="outer",
               
               tags$head(
                 # Include our custom CSS
                 includeCSS("styles.css")#,
                 #includeScript("gomap.js")
               ),
               
               leafletOutput("mymap", width="100%", height="100%"),
               
               # Shiny versions prior to 0.11 should use class="modal" instead.
               absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                             draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                             width = 330, height = "auto",
                             
                             p(),
                             h2('Intervalo de Tiempo'),
                             dateInput("date", "Inicio"),
                             dateInput("date2","Final"),
                             # sliderInput("animation", "Animacion:",0, 5, 1, step = 1, 
                             #             animate=animationOptions(interval=500, loop=T)),
                             selectInput("dataset1", "Precipitacion", c("RT3B42", "GPM")),
                             actionButton('des','descarga'),
                             actionButton('ppt','plot'),
                             actionButton('merge','Generar Pisco del Dia'),
                             actionButton('extraccion','Tablas-Raster'),
                             h4('Realizado por Jose A. Zevallos Ruiz y Cesar Aybar'),
                             h4('contacto: ingjosezr@gmail.com , cesaraybar1994@gmail.com')
                             
               )
           )#,
           
  ))

