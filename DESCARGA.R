

#instalacion de librerias necesarias 
if (require(sp)==F) install.packages('sp',repos='http://cran.us.r-project.org')
if (require(shiny)==F) install.packages('shiny',repos='http://cran.us.r-project.org')
if (require(zoo)==F) install.packages('zoo',repos='http://cran.us.r-project.org')
if (require(rgeos)==F) install.packages('rgeos',repos='http://cran.us.r-project.org')
if (require(gstat)==F) install.packages('gstat',repos='http://cran.us.r-project.org')
if (require(automap)==F) install.packages('automap',repos='http://cran.us.r-project.org')
if (require(rgdal)==F) install.packages('rgdal',repos='http://cran.us.r-project.org')
if (require(hydroGOF)==F) install.packages('hydroGOF',repos='http://cran.us.r-project.org')
if (require(raster)==F) install.packages('raster',repos='http://cran.us.r-project.org')
if (require(gdalUtils)==F) install.packages('gdalUtils',repos='http://cran.us.r-project.org')
if (require(RCurl)==F) install.packages('RCurl',repos='http://cran.us.r-project.org')
if (require(XML)==F) install.packages('XML',repos='http://cran.us.r-project.org')
if (require(RSelenium)==F) install.packages('RSelenium',repos='http://cran.us.r-project.org')
if (require(doSNOW)==F) install.packages('doSNOW',repos='http://cran.us.r-project.org')
if (require(foreach)==F) install.packages('foreach',repos='http://cran.us.r-project.org')
if (require(leaflet)==F) install.packages('leaflet',repos='http://cran.us.r-project.org')
if (require(tcltk)==F) install.packages('tcltk',repos='http://cran.us.r-project.org')


library(shiny)


path0='C:/Control'

runApp(appDir = path0,launch.browser = T)
