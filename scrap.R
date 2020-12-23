w.pm <- w[wRepl.index$w.repl == i]
  PlotField2(field = w.pm, 
             mesh = meshPois, 
             xlim = range(meshPois$loc[,1]), 
             ylim = range(meshPois$loc[,2]),
             zlim = range(-3.5,3.5), 
             MyMain = MyTitle)
  points(x = Loc[MA3$year==Years[i],1],
         y = Loc[MA3$year==Years[i],2], 
         cex = 0.5, 
         col = "black", 
         pch = 16)


library(sp)
library(albersusa)

spdf <- rmapshaper::ms_simplify(usa_sf(), keep = 0.1)
pal <- colorNumeric("Blues", domain = spdf$pop_2014)
epsg2163 <- leafletCRS(
  crsClass = "L.Proj.CRS",
  code = "EPSG:2163",
  proj4def = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs",
  resolutions = 2^(16:7))

leaflet(spdf, options = leafletOptions(crs = epsg2163)) %>%
  addPolygons(weight = 1, color = "#444444", opacity = 1,
    fillColor = ~pal(pop_2014), fillOpacity = 0.7, smoothFactor = 0.5,
    label = ~paste(name, pop_2014),
    labelOptions = labelOptions(direction = "auto"))


  setView(11.965053, 57.70451, 16) %>%
  addTiles() %>%
  addMarkers(11.965053, 57.70451)
leaflet() %>%
crs <- leafletCRS(proj4def = "+proj=utm +zone=35 +south ellps=WGS84 +datum=WGS84")

library(raster)

r <- raster("nc/oisst-sst.nc")
pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(r),
  na.color = "transparent")

leaflet() %>% addTiles() %>%
  addRasterImage(r, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r),
    title = "Surface temp")