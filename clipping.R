library(maptools)  ## For wrld_simpl
library(raster)

## Example SpatialPolygonsDataFrame
data(wrld_simpl)
SPDF <- subset(wrld_simpl, NAME=="Brazil")
class(SPDF)

## Example RasterLayer
r <- raster(nrow=1e3, ncol=1e3, crs=proj4string(SPDF))
plot(r)
r[] <- 1:length(r)

## crop and mask
r2 <- crop(r, extent(SPDF))
plot(r2)
r3 <- mask(r2, SPDF)
plot(r3)
## Check that it worked
plot(r3)
plot(SPDF, add=TRUE, lwd=2)



##############
library("faraway")
data(penicillin)

View(penicillin)
Z <- as(model.matrix(~ 0 + blend, data = penicillin), "Matrix")
dim(Z)
penicillin$ID <- 1:nrow(penicillin)

library("lme4")
data(sleepstudy)
head(sleepstudy)
str(sleepstudy)

sleepstudy$Reaction <- sleepstudy$Reaction / 1000
