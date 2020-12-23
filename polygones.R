rm(list=ls(all=TRUE))
library(sp)
library(rgeos)
library(leaflet) #, lib="C:/Users/eliecer.diaz/Documenöts/R/win-library/3.6")
#install.packages("leaflet", lib="C:/Users/eliecer.diaz/Documents/R/win-library/3.6")

Sr1 = Polygon(cbind(c(2, 4, 4, 1, 2), c(2, 3, 5, 4, 2)))
Sr2 = Polygon(cbind(c(5, 4, 2, 5), c(2, 3, 2, 2)))
Sr3 = Polygon(cbind(c(4, 4, 5, 10, 4), c(5, 3, 2, 5, 5)))
Sr4 = Polygon(cbind(c(5, 6, 6, 5, 5), c(4, 4, 3, 3, 4)), hole = TRUE)
Srs1 = Polygons(list(Sr1), "s1")
Srs2 = Polygons(list(Sr2), "s2")
Srs3 = Polygons(list(Sr4, Sr3), "s3/4")
SpP = SpatialPolygons(list(Srs1, Srs2, Srs3), 1:3)

leaflet(height = "300px") %>% addPolygons(data = SpP)


###########################
# https://towardsdatascience.com/animating-your-data-visualizations-like-a-boss-using-r-f94ae20843e3
###########################
library(gapminder)
library(ggplot2)
install.packages("devtools")
devtools::install_github('thomasp85/gganimate')

transition_time <- function(time, range = NULL) {
  time_quo <- enquo(time)
  require_quo(time_quo, 'time')
  ggproto(NULL, TransitionTime,
          params = list(
            time_quo = time_quo,
            range = range
          )
  )
}


class(gapminder)
dim(gapminder)
colnames(gapminder)
head(gapminder)
ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'Salary', y = 'Visits') +
  transition_time(year) +
  ease_aes('linear')


# Old code
ggplot(mtcars) + 
  geom_boxplot(aes(factor(cyl), mpg, frame = gear))

# New code
ggplot(mtcars) + 
  geom_boxplot(aes(factor(cyl), mpg)) + 
  transition_manual(gear)


#########
library(ggplot2)
install.packages("gganimate")
library(gganimate)

ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_boxplot() + 
  # Here comes the gganimate code
  transition_states(
    gear,
    transition_length = 2,
    state_length = 1
  ) +
  enter_fade() + 
  exit_shrink() +
  ease_aes('sine-in-out')