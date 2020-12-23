rm(list=ls(all=TRUE))

library(geoR)
library(ggplot2)

head(gambia)
dim(gambia)

dim(unique(gambia[, c("x", "y")]))

library(dplyr)

d <- group_by(gambia, x, y) %>%
  summarize(
    total = n(),
    positive = sum(pos),
    prev = positive / total
  )
head(d)

library(sp)
library(rgdal)
sps <- SpatialPoints(d[, c("x", "y")],
  proj4string = CRS("+proj=utm +zone=28")
)
spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))

d[, c("long", "lat")] <- coordinates(spst)
head(d)
# MApping prevalence
library(leaflet)
library(viridis)

pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(d) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~long, lat = ~lat, color = ~ pal(prev)) %>%
  addLegend("bottomright",
    pal = pal, values = ~prev,
    title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))


library(raster)
r <- getData(name = "alt", country = "GMB", mask = TRUE)

pal <- colorNumeric("viridis", values(r),
  na.color = "transparent"
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
    pal = pal, values = values(r),
    title = "Altitude"
  ) %>%
  addScaleBar(position = c("bottomleft"))

# adding altitude to coordinates
d$alt <- raster::extract(r, d[, c("long", "lat")])
head(d)

#1- MESH (KM)
D   <- dist(d[, c("long", "lat")])
range(D)
summary(D)

library(INLA)
coo <- cbind(d$long, d$lat)
mesh <- inla.mesh.2d(
  loc = coo, max.edge = c(0.1, 5), #mean/10, maxquantile*2
  cutoff = 0.01 # from minimum distance
)

#2-BUILD SPDE
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

#3-INDEX SET
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#4-PROJECTION MATRIX
A <- inla.spde.make.A(mesh = mesh, loc = coo)

#PREDICTION
dp <- rasterToPoints(r)
dim(dp)

plot(dp)
ra <- aggregate(r, fact = 5, fun = mean)
dp <- rasterToPoints(ra)
dim(dp)
plot(ra)

coop <- dp[, c("x", "y")]

Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
dim(Ap)
###################################################
#  Stack with data for estimation and prediction
###################################################
# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$positive, numtrials = d$total), #dependent variable
  A = list(1, A),
  effects = list(data.frame(b0 = 1, altitude = d$alt), #covariates
   s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials = NA), # dependent variable
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, altitude = dp[, 3]),
    s = indexs
  )
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

# Formula
formula <- y ~ 0 + b0 + altitude + f(s, model = spde)

res <- inla(formula,
  family = "binomial", Ntrials = numtrials,
  control.family = list(link = "logit"),
  data = inla.stack.data(stk.full),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stk.full)
  )
)

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

#mapping
pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    lng = coop[, 1], lat = coop[, 2], #coordinates
    color = pal(prev_mean) # Z dimesnioion
  ) %>%
  addLegend("bottomright",
    pal = pal, values = prev_mean, # DATA_legend
    title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))

  # rasterize and plot
  r_prev_mean <- rasterize(
  x = coop, y = ra, field = prev_mean,
  fun = mean
)

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
    pal = pal,
    values = values(r_prev_mean), title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))