# Let's read the jeoJson file that is stored on the web with the geojsonio library:
library(geojsonio)
library(mapproj)
spdf <- geojson_read("https://raw.githubusercontent.com/gregoiredavid/france-geojson/master/communes.geojson",  what = "sp")
 
# We can plot that!
library(sp)
plot(spdf, lwd=0.1)
class(spdf)

# head(spdf)
library(broom)
spdf_fortified <- tidy(spdf, region = "code")
 
# Now I can plot this shape easily as described before:
ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group)) +
  theme_void() +
  coord_map()