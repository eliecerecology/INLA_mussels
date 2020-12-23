rm(list=ls(all=TRUE))
########################
#map of the event
########################

library(GISTools)
library(rgdal)
library(raster)
#library(tidyverse)
library(osmdata)
library(OpenStreetMap)
library(sf)
library(ggmap)

mad_map <- openmap(c(60.19121, 24.91830),
                   c(60.18778, 24.92828),
                   type = "osm",
                   zoom = NULL) #,

mad_mapP <- openproj(mad_map, projection="+proj=longlat +ellps=WGS84 +datum=WGS84  ")

plot(mad_mapP)
axis(1)
axis(2)
mtext(side= 1, line= 3, text = "longitude", cex =1.2, font =2)
mtext(side= 2, line= 2, text = "latitude", cex =1.2, font =2)

p <- cbind(runif(40,  24.921, 24.923), runif(40, 60.1889, 60.1897))
points(p)

# NOT RUN {
# data.frame
A = runif(100,0, 0.003)
data = data.frame(x= (24.922 + A/2 + runif(10, 0, 0.0007)),
                  y = (60.188 + A/2 + runif(10, 0, 0.0008)))

#coordinates(data) <- ~x+y

data.xy = data[c("x", "y")]
coordinates(data.xy) <- ~x+y
class(data.xy)
points(data.xy,  cex = 1.5, pch = 21, col="red")

# install.packages("clusteringdatasets")
# blobs <- make_blobs()
# plot(blobs$samples, col=rainbow(3)[blobs$labels])

# library(clusterSim)
# means <- matrix(c(0,7,0,7),2,2)
# cov <- matrix(c(1,0,0,1),2,2)
# grnd <- cluster.Gen(numObjects=60,means=means,cov=cov,model=2,
# numOutliers=8)
# colornames <- c("red","blue","green")
# grnd$clusters[grnd$clusters==0]<-length(colornames)
# plot(grnd$data,col=colornames[grnd$clusters],ask=TRUE)

library(MASS)
# Simulate bivariate normal data
mu <- c(22.923, 62.190)                         # Mean
Sigma <- matrix(c(1, 0.0000000000000000000000001, 0.0000000000000000000000001, 1), 2)  # Covariance matrix

# > Sigma
# [,1] [,2]
# [1,]  1.0  0.1
# [2,]  0.1  1.0
 
# Generate sample from N(mu, Sigma)
bivn <- mvrnorm(500, mu = mu, Sigma = Sigma )  # from Mass package
head(bivn) 
class(bivn)

# Calculate kernel density estimate
bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)   # from MASS package
image(bivn.kde)

bivn <- data.frame(bivn)
head(bivn)
colnames(bivn) <- c("x", "y")
datass.xy = bivn[c("x", "y")]

coordinates(datass.xy) <- ~x + y

class(datass.xy)
datass.xy <- openproj(datass.xy, projection="+proj=longlat +ellps=WGS84 +datum=WGS84  ")

plot(mad_map)

plot(datass.xy,  cex = 1.5, pch = 21, col="red")
axis(1)
axis(2)


######################


# Draw from multi-t distribution without truncation
library (tmvtnorm)
Sigma <- matrix(c(1, .1, .1, 1), 2)  # Covariance matrix
X1 <- rtmvt(n=1000, mean=rep(0, 2), sigma = Sigma, df=2) # from tmvtnorm package
 
t.kde <- kde2d(X1[,1], X1[,2], n = 50)   # from MASS package
col2 <- heat.colors(length(bivn.kde$z))[rank(bivn.kde$z)]
persp3d(x=t.kde, col = col2)

#######################
cbind(rnorm(40) * 2 + 13, rnorm(40) + 48)