rm(list=ls(all=TRUE))

#inla.upgrade(testing = FALSE)
library(ggplot2)
library(ggmap)
library(INLA)
library(sp)
library(dplyr)
library(gstat)
library(rgdal)
library(raster)
library(dismo)
library(splancs)
library(reshape)
library(lattice)

library(inlatools)
library('maps')
library('maptools')
library('mapdata')
data("worldHiresMapEnv")
source("/home/elvi/Documents/MegaSync/RCourseInla/AllRCode/HighstatLibV11.R")
source("HighstatLibV10.R")
load(file = "/home/elvi/Documents/MegaSync/Recruitment/MA1.Rda")

LongLatToUTM<-function(x,y){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS("+proj=utm +zone=35"))
  return(as.data.frame(res))
}

xy <- LongLatToUTM(x = MA1$Lon, y = MA1$Lat)
MA1$Xutm <- xy[,2]
MA1$Yutm <- xy[,3]

#convert UTM to long.lat

MA1$Xkm <- MA1$Xutm / 1000
MA1$Ykm <- MA1$Yutm / 1000
#MA1$Xkm <- MA1$Xutm 
#MA1$Ykm <- MA1$Yutm 
###############
# MAPA clipped
#################
library(rgdal)
DSN <- "/home/elvi/Downloads/costline2/terrain.shp"
ShapeF.utm <- readOGR(dsn = DSN, layer = "terrain")
crs(ShapeF.utm)
#crs(ShapeF.utm) <- "+proj=utm +zone=35N"
plot(ShapeF.utm)
# Sweet!

##########
# Border
##########
imahe  = rgdal::readOGR("/home/elvi/Downloads/costline2/border.shp")
spplot(imahe, main = "Tvärminne region, Finland")
crs(imahe)

#################3
# Exploration
###################
MyVar <- c("nine", "four","two", "one", "zerofive", 
             "Biomass", "sal0", "temp0",
             "recr_aut", "Ykm")

#Mydotplot(MA1[,MyVar])
# REMOVE OUTLIERS
#MA1 <- MA1[!MA1$two > 1500,]
#Pairs
#Mypairs(MA1[,MyVar])
# corvif(MA1[,MyVar])
# MyVar <- c("adul_juv", "Biomass", "sal0", "temp0")
# corvif(MA1[,MyVar])

# Environmental variables
MA1$sal0.std <- MyStd(MA1$sal0)
MA1$temp0.std <- MyStd(MA1$temp0)
MA1$chl.std <- MyStd(as.numeric(MA1$chlsumer))

# Density depdendent
MA1$nine.std <- MyStd(MA1$nine)
MA1$four.std <- MyStd(MA1$four)
MA1$two.std <- MyStd(MA1$two)
MA1$one.std <- MyStd(MA1$one)
MA1$zerofive.std <- MyStd(MA1$zerofive)
MA1$Xstd <- MyStd(MA1$Xkm)
MA1$Ystd <- MyStd(MA1$Ykm)
MA1$Isaeus.std <- MyStd(MA1$Isaeus)
MA1$recr_aut.std <- MyStd(MA1$recr_aut)

Mycontinuos <- c("sal0.std",
            "temp0.std",
            "chl.std",
            "nine.std",
            "four.std",
            "two.std",
            "one.std",
            "zerofive.std",
            "Xstd",
            "Ystd",
           "recr_aut",
           "Isaeus")

MyVarX <- c("sal0.std",
            "temp0.std",
            "sal0",
            "temp0",
            "Isaeus",
            "chl.std",
            "nine.std",
            "four.std",
            "two.std",
            "one.std",
            "zerofive.std",
            "zerofive",
            "Xstd",
            "Ystd",
            "site",
            "Xutm",
            "Yutm",
            "year",
            "recr_aut",
            "recr_aut.std",
            "Isaeus.std",
            "Lon",
            "Lat")

#corvif(MA1[,Mycontinuos])
MyMultipanel.ggp2(MA1,
                    varx = Mycontinuos,
                    vary = "recr_aut",
                    ylab = "recr_aut",
                    addSmoother = FALSE,
                    addRegressionLine = TRUE,
                    addHorizontalLine = FALSE)

MyMultipanel.ggp2(MA1,
                    varx = Mycontinuos,
                    vary = "zerofive",
                    ylab = "zerofive",
                    addSmoother = TRUE,
                    addRegressionLine = TRUE,
                    addHorizontalLine = FALSE)

################
# REDUCE IT
################ 
MA2 <- MA1[, MyVarX]
MA2$year <- factor(MA2$year)
MA3 <- aggregate(.~ site + year, MA2, mean)
MA3$fournine.std <- MA3$four.std + MA3$nine.std

################
# mesh with border optional
################ 
tv_df <- fortify(imahe)
head(tv_df)

# And we convert the coordinate (which are called lat and long,
# but which are UTM coordinates) to Xkm and Ykm
tv_df$Xkm <- tv_df$long # / 1000
tv_df$Ykm <- tv_df$lat #/ 1000

# Extract the UTM coordinates
CoastLine <- tv_df[,c("Xkm", "Ykm")]
head(CoastLine)

# We can plot these:
xyplot(Ykm ~ Xkm,
       data = CoastLine,
       type = "p",
       aspect = "iso",
       xlab = list("Xkm", cex = 1.0),
       ylab = list("Ykm", cex = 1.0),
       col = 1,
       pch = 16)

# To reverse the locations:
N <- nrow(CoastLine)
ReOrder <- N:1
Loc.Reverse <- CoastLine[ReOrder, c("Xkm", "Ykm")]
plot(Loc.Reverse)
# And here is the mesh with a boundary
mesh2 <- inla.mesh.2d(
          loc.domain = CoastLine,
          max.edge = 1000*1.5,
          boundary = inla.mesh.segment(Loc.Reverse))      
plot(mesh2)
points(MA3$Xutm, MA3$Yutm)

#####################################################
#spatial dependence
# What are distances between sites?
#####################################################
Loc <- cbind(MA3$Xutm, MA3$Yutm)
D   <- dist(Loc)/1000 # to get it in km
# Figure 21.5
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D , 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")

##################
# INLA
##################
# A .NON REPLICATED MESH:

# -1 MESH
RangeGuess <- 10 * 1000 # size triangle inside
MaxEdge <- RangeGuess / 5  # 1000 m # Number triangle inside
ConvHull <- inla.nonconvex.hull(Loc, convex = 1000*10) # outer distance
# The convex option puts the boundary of the innerpart
# closer to the points. Saves some computing time
# during this course. Maybe not use it for a paper.
mesh <- inla.mesh.2d(loc = Loc,
                         boundary = ConvHull, 
                         max.edge = c(1,5) * MaxEdge, # The largest allowed triangle edge length.  One or two values.
                         cutoff = MaxEdge/5) # minimum allowed distance between islands.
mesh$n #276

# In case the use of the border
#meshPois <- mesh2

# Figure 
par(mfrow = c(1,1), mar=c(1, 1, 1, 1))
plot(mesh, asp = 1)
points(Loc, col = 2, pch = 16, cex = 1)

#2- Define the SPDE.
spde <- inla.spde2.pcmatern(meshPois,
                               # alpha = 2, constr = TRUE
                               prior.range = c(10 *1000, 0.01), #not lower than 10km
                               prior.sigma = c(1.6, 0.01)) #SIgmaU

#3- Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spdePois$n.spde)

#########################################
#4- PROJECTION: Define the weighting factors a_ik (also called 
#  the projector matrix).
A  <- inla.spde.make.A(mesh, loc = Loc)


#B. REPLICATED MESH
# Make a stack. 
##################
# Ropes
##################

Repl <- as.numeric(MA3$year) 
table(Repl)
NRepl <- length(unique(Repl))
NRepl #Number of replicated random fields we will get

XYear <- model.matrix(~ year, data = MA3)
N <- nrow(MA3)

# 1-mesh done
meshPois <- inla.mesh.2d(loc = Loc,
                         boundary = ConvHull, 
                         max.edge = c(1,5) * MaxEdge, # The largest allowed triangle edge length.  One or two values.
                         cutoff = MaxEdge/5) # minimum allowed distance between islands.
# 2-SPDE replicated in year
# OPTION 1: spde.PoisRepl <- inla.spde2.matern(meshPois, alpha = 2, constr = TRUE)
# spdePois <- inla.spde2.pcmatern(meshPois,
#                                # alpha = 2, constr = TRUE
#                                prior.range = c(10 *1000, 0.01), #not lower than 10km
#                                prior.sigma = c(1.5, 0.01)) #SIgmaU

spde.PoisRepl <- inla.spde2.pcmatern(meshPois,
                               # alpha = 2, constr = TRUE
                               prior.range = c(10 *1000, 0.01), #not lower than 10km
                               prior.sigma = c(1, 0.01)) #SIgmaU


# 3-Index replicated w
wRepl.index <- inla.spde.make.index('w', 
                                 n.spde = meshPois$n,
                                 n.repl = NRepl)

#4-A Projection matrix
A.PoisRepl <- inla.spde.make.A(meshPois, 
                               loc = Loc, 
                               repl = Repl)
#5-covariates:
CovariatesR <- data.frame(
  Intercept   = rep(1, N),
  year = MA3$year,
  site = MA3$site,
  sal0.std = MA3$sal0.std,
  temp0.std = MA3$temp0.std,
  Isaeus.std = MA3$Isaeus.std)

#6- formula
fPoisR <- y ~ -1 + Intercept +
              #year +
              sal0.std + temp0.std + Isaeus.std +
              f(year, model = "iid") +              
              f(w, model = spde.PoisRepl, replicate = w.repl)

fPoisNO<- y ~ -1 + Intercept +
              #year +
              sal0.std + temp0.std + Isaeus.std +
              f(year, model = "iid") #+ f(site, model = "iid")

fPoislope<- y ~ -1 + Intercept +
              #year +
              temp0.std + Isaeus.std +
              f(year, sal0.std, model = "iid") + f(w, model = spde.PoisRepl, replicate = w.repl)

names(inla.models()$likelihood)

StackPoisReplR <- inla.stack(
  tag = "Fit",
  data = list(y = as.integer(MA3$recr_aut)),  
  A = list(1, A.PoisRepl),                  
  effects = list(   
    Covariates = CovariatesR,
    w          = wRepl.index))

prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.1,
                    mean = list(sal0.std = 0, prior = "lognormal",
                                temp0.std = 0, prior = "normal",
                                Isaeus.std = 0, prior = "lognormal"),
                    prec = list(sal0.std = 1,
                                temp0.std = 1,
                                Isaeus.std = 1)
                    )                    

G3R <- inla(fPoisR, 
           family = "poisson", 
           data=inla.stack.data(StackPoisReplR),
           control.fixed = prior.fixed,
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackPoisReplR)))

# G3RNO <- inla(fPoisNO, 
#            family = "poisson", 
#            data=inla.stack.data(StackPoisReplR),
#            control.fixed = prior.fixed,
#            control.compute = list(dic = TRUE, waic = TRUE),
#            control.predictor = list(A = inla.stack.A(StackPoisReplR)))

G3Rslope <- inla(fPoislope, 
           family = "poisson", 
           data=inla.stack.data(StackPoisReplR),
           control.fixed = prior.fixed,
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackPoisReplR)))
summary(G3R)
plot(simulate_iid(sigma =2))
plot(dispersion_check(G3Rslope))


p   <- nrow(G3R$summary.fixed)
muR <- G3R$summary.fitted.values[1:N, "mean"]
E1  <- (as.integer(MA3$recr_aut) - muR) / sqrt(muR)
sum(E1^2) / (N - p + 2) 
# 0.07

# dic2  <- c(G3R$dic$dic, G3RNO$dic$dic, G3Rslope$dic$dic)
# waic2 <- c(G3R$waic$waic, G3RNO$waic$waic, G3Rslope$waic$waic)
# G3R <- G3Rslope
####################
# PREDICTION
####################
N <- nrow(Covariates_P)

library(plyr)
le = 30
# Covariates <- data.frame(
#   Intercept   = rep(1, N),
#   year = MA3$year,
#   sal0.std = MA3$sal0.std,
#   temp0.std = MA3$temp0.std,
#   Isaeus.std = MA3$Isaeus.std)

# Covariates_P <- data.frame(
#     Intercept   = rep(1, le),
#     sal0.std = seq(min(MA3$sal0.std),  max(MA3$sal0.std), length = le),
#     temp0.std = 0, #seq(min(MA3$temp0.std),  max(MA3$temp0.std), length = le),
#     Isaeus.std  = 0,
#     year = 2008
# )

Covariates_P <- ddply(MA3, .(year), summarize,
    Intercept   = rep(1, le),
    sal0.std = seq(min(MA3$sal0.std),  max(MA3$sal0.std), length = le),
    temp0.std = rep(0, le),   #seq(min(MA3$temp0.std),  max(MA3$temp0.std), length = le),
    Isaeus.std  = rep(0, le)  #seq(min(MA3$Isaeus.std),  max(MA3$Isaeus.std), length = le)
)          
#Covariates_P$year <- NULL

Stack_Pred <- inla.stack(
  tag = "Predict",
  data = list(y = NA),  
  A = list(1, A.PoisRepl),                 
  effects = list(
    Covariates = Covariates_P,
    w          = wRepl.index))

StackPoisReplR <- inla.stack(
  tag = "Fit",
  data = list(y = as.integer(MA3$recr_aut)),  
  A = list(1, A.PoisRepl),                  
  effects = list(   
    Covariates = CovariatesR,
    w          = wRepl.index))

stk.full <- inla.stack(StackPoisReplR, Stack_Pred)

prediction <- inla(fPoislope, 
           family = "poisson", 
           data=inla.stack.data(stk.full),
           control.fixed = prior.fixed,
           control.compute = list( dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk.full)))

index <- inla.stack.index(stk.full, tag = "Predict")$data
res <- prediction$summary.fitted.values[index, c("mean", "0.025quant", "0.975quant")]  #210 by 3

# Log-link prediction
pred_mean <- exp(res[, "mean"])
pred_ll <- exp(res[, "0.025quant"])
pred_ul <- exp(res[, "0.975quant"])

# Double check
# X <- as.matrix(Covariates_P) # ok
# beta <- prediction$summary.fixed[,"mean"] # yes
# wpm <- prediction$summary.random$w$mean # yes
# eta <- exp(X %*% beta + A.PoisRepl %*% as.matrix(wpm)) # yes
# cbind(eta, exp(pred_mean))

marginal_sal <- inla.smarginal(G3R$marginals.fixed$sal0.std)
marginal_sal <- data.frame(marginal_sal)
mar_sal <- ggplot(marginal_sal, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[1], "salinity")), y = "Density",
        title = c("Posterior distribution for salinity")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"))
#mar_sal + theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"))

varis <- MA3$sal0
sal <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = varis, y = recr_aut), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = varis, y = pred_mean),
   method="glm", se = FALSE, colour = "black",
   method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = varis, y = pred_ll), 
   method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = varis, y= pred_ul), 
  method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Salinity (psu)",
       y = "Mussel recruit density",
       title = "Mussel recruits in experimental brick collectors",
       subtitle = "Black circles are obeserved values.") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0))

##############
# Temperature
###############
marginal_temp <- inla.smarginal(G3R$marginals.fixed$temp0.std)
marginal_temp <- data.frame(marginal_temp)
mar_temp <- ggplot(marginal_temp, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[2], "temperature")), y = "Density",
  title = c("Posterior distribution for temperature")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"))


vari <- MA3$temp0

temp <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = vari, y = recr_aut), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = vari, y = pred_mean),
   method="glm", se = FALSE, colour = "black",
   method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = vari, y = pred_ll), 
   method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = vari, y= pred_ul), 
  method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Temperature (C)",
       y = "Mussel recruit density",
       title = "Mussel recruits in experimental brick collectors",
       subtitle = "Black circles are obeserved values.") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0))

###########
# ISAEUS
###########
marginal_exp <- inla.smarginal(G3R$marginals.fixed$Isaeus.std)
marginal_exp <- data.frame(marginal_exp)
mar_exp <- ggplot(marginal_exp, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[3], "Wave exposure")), y = "Density",
  title = c("Posterior distribution for wave exposure")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"))


varie <- MA3$Isaeus

we <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = varie, y = recr_aut), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = varie, y = pred_mean),
   method="glm", se = FALSE, colour = "black",
   method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = varie, y = pred_ll), 
   method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = varie, y= pred_ul), 
  method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Wave exposure index (Isaeus)",
       y = "Mussel recruit density",
       title = "Mussel recruits in experimental brick collectors",
       subtitle = "Black circles are obeserved values.") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0))


# insert ggplot code
tiff("/home/elvi/Documents/BayesianAnalysis/figures/Figure1.tiff", units="in", width=15, height=15, res=150)
jpeg("/home/elvi/Documents/BayesianAnalysis/figures/Figure_charts1.jpeg", width=1000, height=1000, res=118.1)
figure <- ggarrange(mar_sal, sal, mar_temp, temp, mar_exp, we,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3)
figure
dev.off()
#######################################
## PLOT################################
library(rgeos)
library(MASS)
library(fields)
library("ggpubr")

PlotField2 <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE, MyMain, ...){
  stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- c(275000, 295000) 
  if (missing(ylim)) ylim <- c(6637799, 6643943) #c(6627799, 6643943)
  
  # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = xlim, 
                              ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  
  # And plot the whole thing
  image.plot(list(x = proj$x, 
                  y = proj$y,
                  z = field.proj), 
             xlim = xlim, 
             ylim = ylim,
             asp = 1,
             add = Add,
             main = MyMain,
             xlab = c("Easting"),
             ylab = c("Northing"),
             cex.lab = 1.5, cex.main = 1.3, cex.axis = 0.8, lwd = 5,
             ...) 
}
# Plot the spatial random field 
# First we are going to create a spatial polygon in
# UTM coordinates for 'land'. We won't explain
# the next block in detail. Just run it.

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤4
Years <- levels(MA3$year)
w     <- G3R$summary.random$w$mean
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤4
#pdf("/home/elvi/Documents/BayesianAnalysis/Project/RcodeINLA/Plot1.pdf")
#tiff("test.tiff", units="in", width=5, height=5, res=300)
# insert ggplot code

jpeg("/home/elvi/Documents/BayesianAnalysis/figures/Figure1_1a.jpeg", width= 800, height=800, res=118.1) #, 
#pdf("/home/elvi/Documents/BayesianAnalysis/figures/Figure1.pdf", width=15, height=15)
#png(file="/home/elvi/Documents/BayesianAnalysis/figures/Figure1_1.png", width=15, height=15)
par(oma=c(10,35,3,10), mar = c(5,5,3,3)) # margin of 4 spaces width at right hand side
set.panel(2,2) # 2X2 matrix of plots

for (i in 1:length(Years)){
  MyTitle <- Years[i]

  w.pm <- w[wRepl.index$w.repl == i]
  PlotField2(field = w.pm, 
             mesh = mesh, 
             xlim = c(280000.8, 300000.8), # 278475.8, 300000.8 
             ylim =  c(6627799, 6643000), #6627799 6643100
             #xlim = c(278475.8, 300000.8), 
             #ylim =  c(6627799, 6643943),
             zlim = range(-3.5,3.5), 
             MyMain = MyTitle)
  points(x = Loc[MA3$year==Years[i],1],
         y = Loc[MA3$year==Years[i],2], 
         cex = 0.8, # size dot
         col = "black", 
         pch = 16) #shape of dot
 plot(ShapeF.utm, add = TRUE, col = "gray")
   
}

dev.off()

###############another graph
# proj <- inla.mesh.projector(meshPois, 
#                               xlim = range(meshPois$loc[,1]), 
#                               ylim = range(meshPois$loc[,2]), 
#                               dims = c(300, 300))

# mean_s <- inla.mesh.project(proj, w.pm)
# df <- expand.grid(x = proj$x, y = proj$y)
# df$mean_s <- as.vector(mean_s)
# ########TRANSBACK#
# mxy <- data.frame(x = df$x, y = df$y)
# mapdata <- mxy
# colnames(mapdata) <- c("x","y")
# coordinates(mapdata) <- ~x+y #similar to SpatialPoints
# proj4string(mapdata) <- CRS("+proj=utm +zone=35") #assign projection and coordinate reference system
# longlats <- spTransform(mapdata, CRS("+proj=longlat")) #transform
# longlats <- as.data.frame(longlats)
# longlats$mean_s <- as.vector(mean_s)
# head(longlats)
# head(df)

# library(viridis)
# library(cowplot)
# gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
#   #geom_tile() + scale_alpha(range = c(-10, 0.5)) +
#   geom_raster() +
#   scale_fill_viridis(na.value = "transparent") +
#   coord_fixed(ratio = 1) + theme_bw()

# set.panel(1,1)

# plot_grid(gmean)
# plot(ShapeF.utm, add= TRUE, col = "grey")
##################
# Ropes
##################

Repl <- as.numeric(MA3$year) 
table(Repl)
NRepl <- length(unique(Repl))
NRepl #Number of replicated random fields we will get

XYear <- model.matrix(~ year, data = MA3)
N <- nrow(MA3)

# 1-mesh done
meshPois <- inla.mesh.2d(loc = Loc,
                         boundary = ConvHull, 
                         max.edge = c(1,5) * MaxEdge, # The largest allowed triangle edge length.  One or two values.
                         cutoff = MaxEdge/5) # minimum allowed distance between islands.
# 2-SPDE replicated in year
# OPTION 1: spde.PoisRepl <- inla.spde2.matern(meshPois, alpha = 2, constr = TRUE)
# spdePois <- inla.spde2.pcmatern(meshPois,
#                                # alpha = 2, constr = TRUE
#                                prior.range = c(10 *1000, 0.01), #not lower than 10km
#                                prior.sigma = c(1.5, 0.01)) #SIgmaU

spde.PoisRepl <- inla.spde2.pcmatern(meshPois,
                               # alpha = 2, constr = TRUE
                               prior.range = c(10 *1000, 0.01), #not lower than 10km
                               prior.sigma = c(1, 0.01)) #SIgmaU


# 3-Index replicated w
wRepl.index <- inla.spde.make.index('w', 
                                 n.spde = meshPois$n,
                                 n.repl = NRepl)

#4-A Projection matrix
A.PoisRepl <- inla.spde.make.A(meshPois, 
                               loc = Loc, 
                               repl = Repl)

#######################################
# 2 ANALYSES benthic zero five: DD
###################################
CovariatesZ <- data.frame(
  Intercept   = rep(1, N),
  year = MA3$year,
  site = MA3$site,
  fournine.std = MA3$fournine.std,
  two.std = MA3$two.std,
  one.std = MA3$one.std,
  recr_aut.std = MA3$recr_aut.std)

StackPoisRepl_zero <- inla.stack(
  tag = "Fit",
  data = list(y = as.integer(MA3$zerofive)),  
  A = list(1, A.PoisRepl),                  
  effects = list(   
    Covariates = CovariatesZ,
    w          = wRepl.index))

fPoisZ <- y ~ -1 + Intercept + 
              fournine.std +
              two.std +
              one.std +
              #recr_aut.std +
              #f(year, model= "iid" ) +
              f(w, model = spde.PoisRepl, replicate = w.repl) 
fPoisZ <- y ~ -1 + Intercept +
              f(year, site, model = "iid") + f(w, model = spde.PoisRepl, replicate = w.repl)
fPoisZ_i <- y ~ -1 + Intercept + 
              #fournine.std +
              #two.std +
              #one.std +
              #recr_aut.std +
              f(year, model= "iid" ) +
              f(w, model = spde.PoisRepl, replicate = w.repl) 
 
prior.fixedZ <- list(mean.intercept = 0, prec.intercept = 0.1,
                    mean = list(fournine.std = 0, prior = "normal",
                                two.std = 0, prior = "normal",
                                one.std = 0, prior = "normal"),
                    prec = list(fournine.std = 10,
                                two.std = 10,
                                one.std = 10)
                    )                    

G3Z <- inla(fPoisZ, #fGam3, >>>>>>>>>>>>>>>>>>>>><FOR MATS PAPER)
           family = "poisson", 
           data=inla.stack.data(StackPoisRepl_zero),
           control.fixed = prior.fixedZ,
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackPoisRepl_zero)))


 dic2  <- c(G3Z$dic$dic, G3Z1$dic$dic)
# waic2 <- c(G3Z$waic$waic, G3Z_i$waic$waic)
#G3Z_i <- G3Z
#summary(G3Z)
p   <- nrow(G3Z$summary.fixed)
muZ <- G3Z$summary.fitted.values[1:N, "mean"]
E1  <- (MA3$zerofive - muZ) / sqrt(muZ)
sum(E1^2) / (N - p)
#0.08 ~ 0.1

Years <- levels(MA3$year)  # FOR PLOTTING
w     <- G3Z$summary.random$w$mean


tiff("/home/elvi/Documents/BayesianAnalysis/figures/FigureV7.tiff", units = "in", width= 17, height=17, res = 300)
par(oma=c(10,1,1,1), mar = c(8,6,8,6)) # 10,35,3,10 margin of 4 spaces width at right hand side
set.panel(2,2) # 2X2 matrix of plots


for (i in 1:length(Years)){
  MyTitle <- Years[i] 

  w.pm <- w[wRepl.index$w.repl == i]
  PlotField2(field = w.pm, 
             mesh = mesh, 
             xlim = c(280000.8, 300000.8), # 278475.8, 300000.8 
             ylim =  c(6627799, 6643000), #6627799 6643100
             #xlim = c(278475.8, 300000.8), 
             #ylim =  c(6627799, 6643943),
             zlim = range(-3.5,3.5), 
             MyMain = NULL)
  points(x = Loc[MA3$year==Years[i],1],
         y = Loc[MA3$year==Years[i],2], 
         cex = 0.8, # size dot
         col = "black", 
         pch = 16) #shape of dot
 plot(ShapeF.utm, add = TRUE, col = "gray")
   
}
dev.off()


##########################3
# New Plot
####################
# PREDICTION
####################
library(plyr)
le = 30
# Covariates <- data.frame(
#   Intercept   = rep(1, N),
#   year = MA3$year,
#   sal0.std = MA3$sal0.std,
#   temp0.std = MA3$temp0.std,
#   Isaeus.std = MA3$Isaeus.std)

# Covariates_P <- data.frame(
#     Intercept   = rep(1, le),
#     sal0.std = seq(min(MA3$sal0.std),  max(MA3$sal0.std), length = le),
#     temp0.std = 0, #seq(min(MA3$temp0.std),  max(MA3$temp0.std), length = le),
#     Isaeus.std  = 0,
#     year = 2008
# )

Covariates_PZ <- ddply(MA3, .(year), summarize,
    Intercept   = rep(1, le),
    fournine.std = seq(min(MA3$fournine.std),  max(MA3$fournine.std), length = le),
    two.std = seq(min(MA3$two.std),  max(MA3$two.std), length = le),
    one.std = seq(min(MA3$one.std),  max(MA3$one.std), length = le),
    recr_aut.std = seq(min(MA3$recr_aut.std),  max(MA3$recr_aut.std), length = le)
)    
#Covariates_P$year <- NULL

Stack_PredZ <- inla.stack(
  tag = "Predict",
  data = list(y = NA),  
  A = list(1, A.PoisRepl),                 
  effects = list(
    Covariates = Covariates_PZ,
    w          = wRepl.index))

StackPoisReplZ <- inla.stack(
  tag = "Fit",
  data = list(y = as.integer(MA3$zerofive)),  
  A = list(1, A.PoisRepl),                  
  effects = list(   
    Covariates = CovariatesZ,
    w          = wRepl.index))

stk.full <- inla.stack(StackPoisReplZ, Stack_PredZ)

prediction <- inla(fPoisZ, 
           family = "poisson", 
           data=inla.stack.data(stk.full),
           control.fixed = prior.fixed,
           control.compute = list( dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(stk.full)))

index <- inla.stack.index(stk.full, tag = "Predict")$data
res <- prediction$summary.fitted.values[index, c("mean", "0.025quant", "0.975quant")]  #210 by 3

# Log-link prediction
pred_mean <- exp(res[, "mean"])
pred_ll <- exp(res[, "0.025quant"])
pred_ul <- exp(res[, "0.975quant"])

# Double check
# X <- as.matrix(Covariates_P) # ok
# beta <- prediction$summary.fixed[,"mean"] # yes
# wpm <- prediction$summary.random$w$mean # yes
# eta <- exp(X %*% beta + A.PoisRepl %*% as.matrix(wpm)) # yes
# cbind(eta, exp(pred_mean))

marg_49 <- inla.smarginal(G3Z$marginals.fixed$fournine.std)
marg_49 <- data.frame(marg_49)
mar_49 <- ggplot(marg_49, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[1], "Adult sizes > 4 mm (density)")), y = "Density",
  title = c("Posterior distribution for adult mussels")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw()


varif <- MA3$fournine
fournine <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = varif, y = zerofive), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = varif, y = pred_mean),
   method="glm", se = FALSE, colour = "black",
   method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = varif, y = pred_ll), 
   method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = varif, y= pred_ul), 
  method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Adult density",
       y = "Mussel recruits in benthos",
       title = "Mussel recruits found in benthos",
       subtitle = "Black circles are obeserved values") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0))

##############
# 2 mm
###############
marginal_two <- inla.smarginal(G3Z$marginals.fixed$two.std)
marginal_two <- data.frame(marginal_two)
mar_two <- ggplot(marginal_two, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[2], "2 mm")), y = "Density",
  title = c("Posterior distribution for mussels of 2 mm")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw()


var2 <- MA3$two

two <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = var2, y = zerofive), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = var2, y = pred_mean),
   method="lm", se = FALSE, colour = "black",
   #method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = var2, y = pred_ll), 
   method="lm", se = FALSE, colour = "red", 
   #method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = var2, y= pred_ul), 
  method="lm", se = FALSE, colour = "red", 
   #method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Mussel density of 2 mm",
       y = "Mussel recruits in benthos",
       title = "Mussel recruits found in benthos",
       subtitle = "Black circles are obeserved values") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0))

###########
# 1 mm
###########
marginal_1 <- inla.smarginal(G3Z$marginals.fixed$one.std)
marginal_1 <- data.frame(marginal_1)
mar_1 <- ggplot(marginal_1, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[3], "1 mm")), y = "Density",
  title = c("Posterior distribution for mussels of 1 mm")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw()


vari1 <- MA3$one

one <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = vari1, y = zerofive), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = vari1, y = pred_mean),
   method="lm", se = FALSE, colour = "black",
   #method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = vari1, y = pred_ll), 
   method="lm", se = FALSE, colour = "red", 
   #method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = vari1, y= pred_ul), 
  method="lm", se = FALSE, colour = "red", 
   #method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Mussel density of 1 mm",
       y = "Mussel recruits in benthos",
       title = "Mussel recruits found in benthos",
       subtitle = "Black circles are obeserved values") +
  #scale_x_continuous(expand = c(5.25, 6.25)) +
  #facet_wrap(~year) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0))

###########
# Recruit ropes
###########

marginal_rr <- inla.smarginal(G3Z$marginals.fixed$recr_aut)
marginal_rr <- data.frame(marginal_rr)
mar_rr <- ggplot(marginal_rr, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(paste(beta[4], "recruit mussels in ropes")), y = "Density",
  title = c("Posterior distribution for mussels in ropes")) +
  geom_vline(xintercept = 0, col = "black") + theme_bw()


variRR <- MA3$recr_aut

rr <- MA3 %>% bind_cols(as_tibble(res)) %>%
  ggplot() +
  geom_point(aes(x = variRR, y = zerofive ), shape = 1, size = 2, stroke = 1.05) +
  
  geom_smooth(aes(x = variRR, y = pred_mean),
   method="glm", se = FALSE, colour = "black",
   method.args = list(family = "poisson")
   ) +
  
  geom_smooth(aes(x = variRR, y = (pred_ll)), 
   method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +

  geom_smooth(aes(x = variRR, y= pred_ul), 
  method="glm", se = FALSE, colour = "red", 
   method.args = list(family = "poisson"),
   linetype = "dashed") +
  
  labs(x = "Mussel density in the experimental brick collector",
       y = "Mussel recruits in benthos",
       title = "Mussel recruits found in experimental brick collector",
       subtitle = "Black circles are obeserved values") +
  #facet_wrap(~year) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0))

tiff("/home/elvi/Documents/BayesianAnalysis/figures/Figure3.tiff", units="in", width=15, height=15, res=150)

figure2 <- ggarrange(mar_49, fournine, mar_two, two, mar_1, one,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3)
figure2
dev.off()


########################
# Hyperparameters
########################
SpFi.w <- inla.spde2.result(inla = G3Z,
                            name = "w",
                            spde = spde,
                            do.transfer = TRUE)

names(SpFi.w)

Kappa <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.kappa[[1]] )

sigmau <- inla.emarginal(function(x) sqrt(x), 
                        SpFi.w$marginals.variance.nominal[[1]] )

r <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.range.nominal[[1]] )

Kappa
sigmau
r