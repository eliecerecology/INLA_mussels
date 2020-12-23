

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
source(file = "HighstatLibV10.R")
#setwd("C:/Users/eliecer.diaz/Documents/Recruitment/")# CGI Windows7 work
#source(file = "C:/Users/eliecer.diaz/Documents/Recruitment/Rcode/HighstatLibV9.R")

MA <- read.csv("/home/elvi/Documents/MegaSync/Recruitment/FinalmenteData.csv", 
                   header = TRUE,
                   sep = ",")
##############
#Houskeeping
##############
MA$year <- factor(MA$year)
MA1 <- MA %>%
  mutate(site = recode(site, 
     'Hyljekivi' = "Hyljekivi",
     'Hyljekivi/Yttergrund' = "Hyljekivi",
     'Hyljekivi/yttergrund' = "Hyljekivi")
    )

MA1$Lat <- as.character(MA1$Lat)
MA1$Lon <- as.character(MA1$Lon)

#########################
#Replace the comma
#########################
MA1$Lat <- (gsub(",", ".", gsub("\\.", "", MA1$Lat)))
MA1$Lat <- (gsub("°", "", MA1$Lat))
MA1$Lat <- as.numeric(MA1$Lat)

MA1$Lon <- (gsub(",", ".", gsub("\\.", "", MA1$Lon)))
MA1$Lon <- (gsub("°", "", MA1$Lon))
MA1$Lon <- as.numeric(MA1$Lon)

MA1$X <- NULL
MA1$X.1 <- NULL

####################
# 1. outliers
####################
MyVar <- c("recr_aut")
#MyDotplot.ggp2(Z = MA1, varx = MyVar)
MA1 <- MA1[!MA1$zerofive > 3000,]


##############################
# 2. how many obs per site
##############################
# table(MA1$site) # 20 per site
# table(MA1$site, MA1$year)
# table(MA1$year, MA1$chlapr)
# table(MA1$year, MA1$chlsumer)

##########
MA1$chlapr <- factor(MA1$chlapr) 
MA1$chlsumer <- factor(MA1$chlsumer) 

##########################
# Salinity as a factor
#########################
library(plyr)
aver_dal <- ddply(MA1, .(year), summarize, mean=mean(sal0))
MA2 <- merge(aver_dal, MA1, by ="year")

MA2$mean <- factor(MA2$mean)

MA1 <- MA2
#######################################

# ##############################
# # Secchi as factor
# #############################3
aver_secchi <- ddply(MA1, .(site), summarize, mean_secchi = mean(secchi))

MA2 <- merge(aver_secchi, MA1, by ="site")

MA1 <- MA2


###################################################################
# prepare UTM coordinates matrix
###################################################################

utmcoor <- SpatialPoints(cbind(MA1$Lon,MA1$Lat), 
                         proj4string = CRS("+proj=longlat"))
# Now we have UTM Easting and Northing.
# COnvert these guys to latitude/longitude
longlatcoor <- spTransform(utmcoor, CRS("+proj=utm +zone=35"))
MA1$Xloc <- longlatcoor$coords.x1
MA1$Yloc <- longlatcoor$coords.x2
# dev.new()
# MA_04 <- MA1[MA1$year == "2004", ]
# MA_05 <- MA1[MA1$year == "2005", ]
# MA_06 <- MA1[MA1$year == "2006", ]
# MA_07 <- MA1[MA1$year == "2007", ]
# M <- c("MA_04", "MA_05", "MA_06", "MA_07"

par(mfrow = c(1, 4), cex.lab = 1.5, mar = c(6,5,2,2))
i = MA_04
xyplot(Yloc ~ Xloc,
            aspect = "iso",
            data = i,
            col = 1,
            cex  = i$zerofive/max(i$zerofive),
            pch = 16)

i = MA_05
xyplot(Yloc ~ Xloc,
            aspect = "iso",
            data = i,
            col = 1,
            cex  = i$zerofive/max(i$zerofive),
            pch = 16)
i = MA_06
xyplot(Yloc ~ Xloc,
            aspect = "iso",
            data = i,
            col = 1,
            cex  = i$zerofive/max(i$zerofive),
            pch = 16)
i = MA_07
xyplot(Yloc ~ Xloc,
            aspect = "iso",
            data = i,
            col = 1,
            cex  = i$zerofive/max(i$zerofive),
            pch = 16)
    


##################################################
# Start of the analysis.

# We will implement the following models.
# Right now we only have the knowledge to fit models
# 1 and 2. We will show a powerpoint presentation
# before applying model 3.

# 1:  There is no spatial-temporal correlation
# 2:  The spatial random field is constant over time
# 3:  The spatial random field is correlated between
#     consecutive years. 
# We have no idea about interactions.

#Model 1: There is no spatial and no spatial-temporal 
#         correlation
#         This is a simple Poisson GLM

# Standardize the covariates
MA1$sal0.std <- MyStd(MA1$sal0)
MA1$temp0.std <- MyStd(MA1$temp0)
MA1$Biom.std <- MyStd(MA1$Biomass)
MA1$four.std <- MyStd(MA1$four)
MA1$two.std <- MyStd(MA1$two)
MA1$Xloc.std <- MyStd(MA1$Xloc)
MA1$Yloc.std <- MyStd(MA1$Yloc)

MyVar <- c("temp0.std",
        "Biom.std",
        "sal0.std",
        #"Isa.std",
        #"four.std",
        #"two.std", 
        #"sec.std", #,
        #"year",
        #"Yloc.std",
        "Xloc.std")

corvif(MA1[,MyVar])    

#  GVIF
# sal0.std   1.822622
# temp10.std 1.919912
# four.std   1.382622
# two.std    1.756164
# sec.std    1.828681
# Xloc.std   1.424581

# g2 <- ggplot(data=MA1, aes(x=sal0, y= zerofive)) +
#     geom_point(stat="identity") + facet_wrap( ~ mean, nrow = 1)
# g2
# pairs(MA1[, MyVar], 
#       lower.panel = panel.cor)

N <- nrow(MA1)
MA1$Xkmr <- MA1$Xloc + rnorm(N, 0, 0.0001)
MA1$Ykmr <- MA1$Yloc + rnorm(N, 0, 0.0001)

# To get an estimiation of Theta parameter
# library(glmmADMB)
# NB <- glmmadmb(zerofive ~ mean*temp0.std* +
#                          Biom.std*mean,
#                data = MA1, family= "nbinom")

# summary(NB)
library(mgcv)
GAM4 <- gamm(zerofive ~  s(temp0.std, by = year, bs = "cr") +
                         Biom.std*year , # by = Region), # bs = "cr")  
                 data = MA1,
                family = poisson,# neg.bin(theta = c(1.9)), #1.64Gamma(link = log),
                correlation = corExp(form = ~ Xkmr + Ykmr, nugget = TRUE))

# MA1$Isaf <- factor(MA1$Isaeus)
GAM5 <- gamm(zerofive ~  s(temp0.std, by = mean, bs = "cr") +
                        s(Biom.std, by= mean, bs = "cr"),
                         data = MA1,
                family = neg.bin(theta = c(2)), #1.64Gamma(link = log),
                correlation = corExp(form = ~ Xkmr + Ykmr, nugget = TRUE))


# save(GAM4, file = "GAM4NB.rda")
# save(GAM5, file = "GAM5NB.rda")
load("GAM4NB.rda")
load("GAM5NB.rda")
#---- Model ---> M4
E2 <- resid(GAM5$gam, type = "pearson")
F2 <- fitted(GAM5$gam)
MA1$F2 <- F2
plot(F2, E2, xlab = "Fitted", ylab = "Residuals")
abline(h = 0, v = 0)

MyVar <- c("mean",
                "temp0.std",
                "Biom.std",
                  "F2")
                     
MA1$E2 <- E2
MyMultipanel.ggp2(Z = MA1, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = T)
library(mgcv)
g1 <- gam(E2 ~ s(Biom.std), data = MA1)
g2 <- gam(E2 ~ s(F2), data = MA1)
g3 <- gam(E2 ~ s(temp0.std), data = MA1)
plot(g3)
summary(g3)

N    <- nrow(MA1)
sum(E2^2) / (N - 2)
#146 --> 142--->117---101---Gam5----0.94---0.95 

##################################
##############PREDICTION
##################################
library(plyr)

# FORMULA= s(sal0.std) + s(temp10.std) + sec.std + year

GAM4 <- gamm(zerofive ~  s(temp0.std, by = mean, bs = "cr") +
                 Biom.std*mean, # by = Region), # bs = "cr")  
                 data = MA1,
                family = neg.bin(theta = c(2.3)), #1.64Gamma(link = log),
                 correlation = corExp(form = ~ Xkmr + Ykmr, nugget = TRUE))

MyData <- ddply(MA1, .(mean), summarize,
                temp0.std = seq(min(temp0.std), max(temp0.std), length = 25),
                Biom.std = seq(min(Biom.std), max(Biom.std), length = 25)
                )

P <- predict(GAM5$gam, newdata = MyData, se = TRUE) #, se = TRUE, type = "response"

MyData$mu <- exp(P$fit )
MyData$ub <- exp(P$fit + 1.96 * P$se.fit )
MyData$lb <- exp(P$fit - 1.96 * P$se.fit )
dev.new()
q <- ggplot()
q <- q + geom_point(data = MA1, 
                    aes(y = zerofive, x = Biom.std),
                    shape = 16, 
                    size = 3)
q<- q + xlab("Standarized (Biom)") + ylab("Recruits (counts)")
q <- q + theme(text = element_text(size=15)) + theme_bw()

q <- q + geom_line(data = MyData, 
                   aes(x = Biom.std, y = mu), 
                   colour = "black")

q <- q + geom_ribbon(data = MyData, 
                     aes(x = Biom.std, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)

q <-q + facet_grid( .~ mean, scales = "fixed")
q

dev.new()
p <- ggplot()
p <- p + geom_point(data = MA1, 
                    aes(y = zerofive, x = temp0.std),
                    shape = 16, 
                    size = 3)
p <- p + xlab("Standarized (TEMP)") + ylab("Recruits (counts)")
p <- p + theme(text = element_text(size=15)) + theme_bw()

p <- p + geom_line(data = MyData, 
                   aes(x = temp0.std, y = mu), 
                   colour = "black")

p <- p + geom_ribbon(data = MyData, 
                     aes(x = temp0.std, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)

p <- p + facet_grid( ~ mean, scales = "fixed")
p

library(gstat)
library(sp)
mydata <- data.frame(E2, 
                    Xkm = MA1$Xloc/1000,
                    Ykm = MA1$Yloc/1000)
coordinates(mydata) <- c("Xkm", "Ykm")

V1 <- variogram(E2 ~ 1, 
                data = mydata, 
                cressie = TRUE)

V2 <- variogram(E2 ~ MA1.Xloc + MA1.Yloc , 
                data = mydata, 
                cressie = TRUE,
                alpha=c(0, 45, 90, 135))



plot(V2,
     type = "o",
     xlab = list(label = "Distance", cex = 1.5),
     ylab = list(label = "Cressie's semivariance/dissimilarity", cex = 1.5),
     col = 1, pch = 16)

plot(V1,
     type = "o",
     xlab = list(label = "Distance", cex = 1.5),
     ylab = list(label = "Cressie's semivariance/dissimilarity", cex = 1.5),
     col = 1, pch = 16)

##############################################
#Model 2: There is spatial correlation and it is constant over time
#         This is a Poisson GLM + spatial correlation.
#         See also the La Palma exercise


# We will apply the following model.

# zerofive_i ~ Poisson(mu_i)
# E(zerofive) = mu_i
# var(zerofive) = mu_i
# mu_i = exp(Intercept + s(sal0.std_i) + s(temp10.std_i) +
                            # s("secchi.std") + u_i

# u_i ~ N(0, SIGMA)
# Use the SPDE approach to estimate SIGMA.

###################################################
# We will implement the following 8 steps.
# 1. Make a mesh.
# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
# 3. Define the SPDE.
# 4. Define the spatial field.
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.
# 6. Specify the model formula in terms of the 
#    response variable, covariates and the 
#    spatial correlated term.
# 7. Run the spatial model in INLA.
# 8. Inspect the results.

########################
#1. Make a mesh.
#   Step 1 of making a mesh:  Get a 
#   sense for the distribution of 
#   distances between sampling locations. 
Loc <- cbind(MA1$Xloc, MA1$Yloc)
#what are the distances between the points?

# library(geosphere)
D <- dist(Loc)
#D <- distm(Loc)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (m)",
     ylab = "Frequency")


plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (m)",
     ylab = "Cumulative proportion")

#########################################################
# We need a grid on top of our sampling points
#########################################################
library(INLA)
# names(inla.models()$likelihood)

bound <- inla.nonconvex.hull(Loc)
mesh <- inla.mesh.2d(boundary = bound, 
                      max.edge = c(1500, 1000), #max triangle edge
                      cutoff   = 1000) #minimum distance allowed between points
mesh$n
# 694
par(mfrow=c(1,1), mar=c(0,0,2,0))
plot(mesh)
points(Loc, col = 1, pch = 16, cex = 1.5)

#################################
# Step 2. Define the weighting factors a_ik (also called 
#         the projector matrix).
A3 <- inla.spde.make.A(mesh, loc = Loc)
dim(A3)  #595 694
############################################
# Step 3. Define the SPDE.
spde <- inla.spde2.matern(mesh, alpha = 2)

w.index <- inla.spde.make.index(
                   name    = 'w', 
                   n.spde  = spde$n.spde,
                   n.group = 1,
                   n.repl  = 1)

# Step 5.	Make a stack. 
# This is a confusing step. We discussed it
# in detail in Chapter 13.

# Make the X matrix
N <- nrow(MA1) #595

X <- data.frame(Intercept= rep(1,N), 
               mean = MA1$mean,
               temp0.std = MA1$temp0.std,
               Biom.std = MA1$Biom.std
               )
# And here is the stack.
StackFit <- inla.stack(
             tag  = "Fit",
             data = list(y = MA1$zerofive),  
	         A    = list(A3, 1),                      
	         effects = list(                 
	              w = w.index,           #Spatial field  
	              X = as.data.frame(X))) #Covariates

#6.	Specify the model formula in terms of the 
#   response variable, covariates and the 
#   spatial correlated term.

# These are the models that we will fit:
# Y_i ~ Poisson(mu_i)
# E(Y_i) = mu_i
# Model 1: log(mu_i) = Covariate stuff
# Model 2: log(mu_i) = Covariate stuff + u_i


inla.setOption(scale.model.default = TRUE)
thet = 3
hyper.prec = list(theta = list(
                          prior = "pc.prec",
                          param = c(thet, 0.01)))

f2 <- y ~ -1 + Intercept + f(temp0.std, mean, scale.model = TRUE, model = "rw1") +
                           + Biom.std*mean +  f(w, model = spde)
f3 <- y ~ -1 + Intercept + f(temp0.std, model = "rw1", constr = FALSE) +
                           + Biom.std + mean + f(w, model = spde)
f4 <- y ~ -1 + Intercept + temp0.std + Biom.std + mean + f(w, model = spde)
#############################################
# 7. Run the spatial model in INLA.
# First we run the model without spatial dependency (I2),
# then the model with spatial dependency (I3).
library(splines)
inla.setOption(num.threads = 4) 

I2 <- inla(f2,
           family = "poisson",
           #hyper = list(theta = list(prior="pc.prec", param=c(1.84,0.01))), 
           data = inla.stack.data(StackFit),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.predictor = list(A = inla.stack.A(StackFit)),
           verbose = FALSE)


I3 <- inla(f3,
           family = "nbinomial", 
           data = inla.stack.data(StackFit),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.predictor = list(A = inla.stack.A(StackFit)),
           verbose = FALSE)

I4 <- inla(f4,
           family = "nbinomial", 
           data = inla.stack.data(StackFit),
           control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
           control.predictor = list(A = inla.stack.A(StackFit)),
           verbose = FALSE)
# First validation
# And compare the two models with DICs and WAICs
dic2  <- c(I2$dic$dic, I3$dic$dic, I4$dic$dic)
waic2 <- c(I2$waic$waic, I3$waic$waic, I4$waic$waic)
Z     <- cbind(dic2, waic2)
rownames(Z) <- c("Smoother_Poisson GLM + SPDE", 
                 "SmootherNB GLM + SPDE",
                  "NB + SPDE")
Z # smoother NF + SPDE

# For Poisson!
Fit1 <- I2$summary.fitted.values[1:595, "mean"]
E1   <- (MA1$zerofive - Fit1) / sqrt(Fit1)
N    <- nrow(MA1)
sum(E1^2) / (N - 3)
# Poisson 39.31! huge!

theta = 5.7
Fit3 <- I3$summary.fitted.values[1:595, "mean"]
theta <- I3$summary.hyperpar[1, "mean"]
E1   <- (MA1$zerofive - Fit3) / sqrt((Fit3 + Fit3^2) / theta)
N    <- nrow(MA1)
sum(E1^2) / (N - 3)
# 0.90 # minimal underdistpersion

Fit4 <- I4$summary.fitted.values[1:595, "mean"]
theta <- I4$summary.hyperpar[1, "mean"]
E1   <- (MA1$zerofive - Fit4) / sqrt((Fit4 + Fit4^2) / theta)
N    <- nrow(MA1)
sum(E1^2) / (N - 3)
# 1.06

MyData <- ddply(MA1, .(mean), summarize,
                temp0.std = seq(min(temp0.std), max(temp0.std), length = 25),
                zerofive = seq(min(zerofive), max(zerofive), length = 25),
                Biom.std = rep(0, 25)
                )

# P <- predict(GAM5$gam, newdata = MyData, se = TRUE) #, se = TRUE, type = "response"

Fit1 <- I3$summary.fitted.values[1:595, "mean"]
MA1$Fitted1 <- Fit1
MA1$Fit1.025 <- I3$summary.fitted.values[1:595, "0.025quant"]
MA1$Fit1.975 <- I3$summary.fitted.values[1:595, "0.975quant"]



q <- ggplot()
q <- q + geom_point(data = MA1, 
                    aes(y = zerofive, x = Biom.std),
                    shape = 16, 
                    size = 3)
q <- q + xlab("Standarized (Biom)") + ylab("Recruits (counts)")
q <- q + theme(text = element_text(size=15)) + theme_bw()

q <- q + geom_line(data = MA1, 
                   aes(x = Biom.std, y = Fitted1), 
                   colour = "black")

q <- q + geom_ribbon(data = MA1, 
                     aes(x = Biom.std, 
                         ymax = Fit1.975, 
                         ymin = Fit1.025),
                     alpha = 0.5)

q <-q + facet_grid( .~ mean, scales = "fixed")
q
dev.new()
q <- ggplot()
q <- q + geom_point(data = MyData, 
                    aes(y = zerofive, x = temp0.std),
                    shape = 16, 
                    size = 3)
q<- q + xlab("Standarized (Temp 0)") + ylab("Recruits (counts)")
q <- q + theme(text = element_text(size=15)) + theme_bw()

q <- q + geom_line(data = MA1, 
                   aes(x = temp0.std, y = Fitted1), 
                   colour = "black")

q <- q + geom_ribbon(data = MA1, 
                     aes(x = temp0.std, 
                         ymax = Fit1.975, 
                         ymin = Fit1.025),
                     alpha = 0.5)

q <-q + facet_grid( .~ mean, scales = "fixed")
q
summary(I3)
save(I3, file = "NB_GLM_SPDE20July.rda")

library(ggregplot)
Efxplot(list(I2, I3))

Beta1 <- I2$summary.fixed[, c("mean",
                              "0.025quant",  
                              "0.975quant")] 
print(Beta1, digits = 3)

Beta3 <- I3$summary.fixed[, c("mean",
                              "0.025quant",  
                              "0.975quant")] 
print(Beta3, digits = 3)


ggField(I3, mesh, Groups = 1) +
  scale_fill_brewer(palette = "Blues") 

dev.new()

ggField(I3, mesh, Groups = 1) +
  scale_fill_brewer(palette = "Blues") 
dev.new()
ggField(I2, mesh, Groups = 1) +
  scale_fill_brewer(palette = "Blues") 



# How does the spatial correlation look like?

####################################################
# We finally present the spatial component, the wks. 
# Their posterior mean values can be obtained via

w.pm <- I3$summary.random$w$mean 
w.sd <- I2$summary.random$w$sd
# This is the spatial field calculated at all the mesh points

# There are various ways to plot the spatial 
# field. One option is to use the function 
# inla.mesh.projector; it creates a lattice 
# using the mesh and specified ranges. 
w.proj <- inla.mesh.projector(mesh) 

# The function inla.mesh.project can then 
# be used to project the posterior mean 
# values on this grid. By default a lattice 
# of 100 by 100 is used.

w.pm100_100 <- inla.mesh.project(w.proj, w.pm)
w.sd100_100 <- inla.mesh.project(w.proj, w.sd)

# This w.pm100_100 is of dimension 100 by 100 
# and is a projection (interpolation and extrapolation) 
# of the random field w. We can use the levelplot 
# function from the lattice package to plot w.pm100_100


grid <- expand.grid(x = w.proj$x*1000, 
                    y = w.proj$y*1000)
grid$z <- as.vector(w.pm100_100)    # RANDOM FIELD           
plot.wpm <- levelplot(z ~ x * y,
          data = grid, 
          scales = list(draw = TRUE),
          xlab = list("Longitude", cex = 1.5),
          ylab = list("Latitude", cex = 1.5),
          main = list("Posterior mean spatial random field", cex = 1.5))

# And do the same for the posterior standard deviation 
grid$z <- as.vector(w.sd100_100)               
plot.wsd <- levelplot(z ~ x * y,
          data = grid, 
          scales = list(draw = TRUE),
          xlab = list("Longitude", cex = 1.5),
          ylab = list("Latitude", cex = 1.5),
          main = list("Posterior sd spatial random field", cex = 1.5))
#And plot both graphs
print(plot.wpm)
names(MA1)

win.graph()  #This is a Mac command. Use  win.graph() for Windows
print(plot.wsd)












# The SDs are slightly large....but at the larger values of
# w.pm (where w.pm +/ 2 * sd does not contain 0) we seem 
# to have a spatial effect.


# And some final info on the parameters
SpFi.w <- inla.spde2.result(inla = I3,
                            name = "w",
                            spde = spde,
                            do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.kappa[[1]] )

sigmau <- inla.emarginal(function(x) sqrt(x), 
                        SpFi.w$marginals.variance.nominal[[1]] )

r <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.range.nominal[[1]] )

###########
THis is the way
#tau <- I1$marginals.hyperpar$`Precision for the Gaussian observations`
#sigma <- inla.emarginal(function(x) (1/sqrt(x)), tau)
#sigma



c(Kappa, sigmau, r)

# MOre 
# Hyperparameters
SpFi.w <- inla.spde2.result(inla = IM2b,
                            name = "w",
                            spde = spde,
                            do.transfer = TRUE)

names(SpFi.w)# HERE THE LIST

Kappa <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.kappa[[1]] )

sigmau <- inla.emarginal(function(x) sqrt(x), 
                        SpFi.w$marginals.variance.nominal[[1]] )

r <- inla.emarginal(function(x) x, 
                        SpFi.w$marginals.range.nominal[[1]] )

Kappa
sigmau
r


##############################
# correlation structure

#Show correlation structure
D     <- as.matrix(dist(mesh5$loc[,1:2]))
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance (km)", 
     ylab = "Correlation",
     xlim = c(0, 200))







###########################
# simulation model
###########################
set.seed(12345)
NSim <- 1000

# Carry out the simulation
Sim <- inla.posterior.sample(n = NSim, result = I3)

# X matrix
Xm  <- model.matrix(~ Biom.std*mean, data = MA1)
N    <- nrow(MA1)

Ysim <- matrix(nrow = N, ncol = NSim)
dim(Ysim) #595 1000

mu.i <- matrix(nrow = N, ncol = NSim)
Xm   <- as.matrix(Xm)
Am   <- as.matrix(A3)
dim(Am) #595 694
dim(Xm) #595 8
mesh$n
dim(A3)
for (i in 1: NSim){
  Betas    <- Sim[[i]]$latent[c(987, 988, 989)]
  wk       <- Sim[[i]]$latent[416:986] #986 -  = 570 vertices
  dim(wk)
  eta      <- Xm %*% Betas + Am %*% wk #181x3 @ 3 x 1 + 181 x 571 @
  mu.i[,i] <- exp(eta)
  Ysim[,i] <- rpois(n = nrow(OCW), lambda = mu.i[,i])
}







####################################################### 
# Spatial-temporal AR1 Poisson model slightly better
#######################################################

Group4 <- as.integer(MA1$year)
table(Group4)

# We have 4 groups
NGroup <- length(unique(Group4))
NGroup

# Make the projector matrix
A3 <- inla.spde.make.A(mesh, 
                       group = Group4, # <----time (extra)
                       loc = Loc)

dim(MA1) # 595  13
mesh$n   # 694
dim(A3)   #595 rows by 2776 (694*4) columns

w3 <- inla.spde.make.index(name = 'w',  #parameter?
                             n.spde = mesh$n, # 694 vertices
                             n.group = NGroup) # 4 years

w3$w.group 
length(w3$w.group)
 694*4

dim(X)
# And make the stack

Stack3 <- inla.stack(
             tag = "Fit",
             data = list(y = MA1$zerofive),  
	         A = list(A3, 1),    #A3 has time              
	         effects = list(                 
	              w3,         
	              X))
#################################
# Autoregressive correlation
#################################

f5 <- y ~ -1 + Intercept +
       f(temp0.std, mean, scale.model = TRUE, model = "rw1") +
        + Biom.std * mean +
        f(w,
        model = spde,
        group = w.group, # YEAR 
        control.group = list(model='ar1'))

I5 <- inla(f5,
           family = "nbinomial", 
           data = inla.stack.data(Stack3),
           control.compute = list(waic = TRUE, dic = TRUE),  
           control.predictor = list(A = inla.stack.A(Stack3)) # spatial field
           )

# And compare the two models with DICs and WAICs
dic2  <- c(I3$dic$dic, I5$dic$dic)
waic2 <- c(I3$waic$waic,I5$waic$waic )
Z     <- cbind(dic2, waic2)
rownames(Z) <- c("Spatial NB GLM",
                 "Spatial-temporal AR1 NB model")
Z # improves Spatial-temporal
w <- I5$summary.random$w$mean #all ws

wproj <- inla.mesh.projector(mesh) 
w4 <- list()
for (i in 1:4){
   w4[[i]] <- inla.mesh.project(
                      wproj, 
                      w[w3$w.group == i])
}

# This is the posterior mean of the spatial field!  
library(lattice)
require(gridExtra)
do.call(function(...) grid.arrange(..., nrow = 3),
 lapply(w4, 
        levelplot, 
        xlab='', 
        ylab='',
        col.regions=topo.colors(16), 
        scale=list(draw=FALSE)))

w <- I5$summary.random$w$mean
A3m <- as.matrix(A3)
v <- A3m %*% w
length(v)
dim(A3m)

# Which we can plot...
dev.new()
MyCex <- 3 * sqrt(abs(v) /max(v))
MyCol <- c(1,2)[as.numeric(v >0) + 1]
MyPch <- c(1,16)[as.numeric(v >0) + 1]
xyplot(Yloc ~ Xloc | factor(year),
       data = MA1,
       layout = c(4,2),
       aspect = "iso",
       pch = MyPch,
       cex = MyCex,
       col = MyCol,
       xlab = list(label = "X-coordinate", cex = 1.5),
       ylab = list(label = "Y-coordinate", cex = 1.5))



# Here is the correlation between the
# observed data and the fitted values.

idat <- inla.stack.index(Stack3, 'Fit')$data
idat #This is an index telling us which
     #rows in the stack correspond to the 
     #data
Fit5 <- I5$summary.linear.predictor$mean[idat]
cor(MA1$zerofive, Fit5)
# 0.77 => Accuracy (somehow)

# A correlation coefficient is not that usefull 
# for a data set with lots of zeros.

# Here is the phi
phi <- I5$summary.hy[5,"mean"] #the temporal autocorrelation value
phi # 0.80

# Here is the posterior distrbution of phi
rf <- inla.spde2.result(I5, 'w', spde, do.transf=TRUE)

plot(I5$marginals.hyper[[5]], 
     type ='l',
     xlab = expression(phi), 
     ylab ='Density')
abline(v = phi, col = 2)













str(MA1)
table(MA1$chlapr)
table(MA1$site)
table(MA1$year)
table(MA1$chlsumer)

####################
# Standarize
####################
MyVar <- c("Isaeus", "nine", "four",
         "two", "one", "adul_juv", 
         "Biomass", "secchi", "sal0",
         "temp0", "temp10")

MA2 <- MA1 %>% mutate_each_(funs(scale(.) %>% as.vector), 
                             vars= MyVar)

#########################################
################### 3D
library(raster)
library(rasterVis)
library(dplyr)
names(MA2)
#########################################
DA1 <- MA2 %>% 
    select(site, year, Lon, Lat, zerofive, sal0, temp0) 

DA01 <- DA1[DA1$year == 2004,]
DA02 <- DA1[DA1$year == 2005,]
DA03 <- DA1[DA1$year == 2006,]
DA04 <- DA1[DA1$year == 2007,]


DA1 <- DA01 %>% 
    select(Lon, Lat, zerofive) 

DA4 <- DA04 %>% 
    select(Lon, Lat, zerofive) 

library(akima)
library(plotly)

s = interp(x = DA03$Lon, y = DA03$Lat, z = DA03$zerofive, duplicate = TRUE)
p1 <- plot_ly(x = s$x, y = s$y, z = s$z) %>% add_surface() %>%
    layout(
    title = "Spatial recruitment visualization",
    scene = list(
      xaxis = list(title = "Lon"),
      yaxis = list(title = "Lat"),
      zaxis = list(title = "Recruits")
    ))

s = interp(x = DA04$Lon, y = DA04$Lat, z = DA04$sal0, duplicate = TRUE)
p3 <- plot_ly(x = s$x, y = s$y, z = s$z) %>% add_surface(type = "surface", colorscale = list(c(0, 1), c("tan", "blue"))) %>%
    layout(
    title = "Spatial recruitment visualization",
    scene = list(
      xaxis = list(title = "Lon"),
      yaxis = list(title = "Lat"),
      zaxis = list(title = "Recruits")
    ))
p3
p <- subplot(p1, p3)
p 

#####################3
# Contour
########################

library(latticeExtra) 

DA04$z <- with(DA04, Lon * Lat) 
# showing data points on the same color scale 
dev.new()

library(rgdal)
library(raster)
#data.shape<-readOGR(dsn="C:/Users/eliecer.diaz/Documents/Recruitment/iho/",layer="FinalMap")
#shp <- shapefile("C:/Users/eliecer.diaz/Documents/Recruitment/iho/FinalMap.shp")
#crs(shp) <- CRS('+init=EPSG:4326')
dev.new()
levelplot(temp0 ~ Lon * Lat, DA04, panel = panel.levelplot.points, cex = 1.2) +
    layer_(panel.2dsmoother(..., n = 200)) 


####################
#####contour 2
####################

plot(DA04$Lon,  DA04$Lat)
z <-  DA04$sal0

#assign lat long coordintes to a variable coor1 using coulmn binding
coor1 <- cbind ( DA04$Lon, DA04$Lat) #x is longitude

dat <- list(X=coor1, Y=z)

thnS <- Tps(dat$X,dat$Y) #thin plate spline function
surface(thnS)
title("Sal0")

###############################################
# contour 3 ---IDW---Inverse distance weight
###############################################
library(gstat)
d = DA04
coordinates(d) <- ~Lon + Lat
summary(d)
Lon.range <- as.numeric( c(23.1372, 23.3031))
Lat.range <- as.numeric( c(59.7575, 59.8640))

grid <- expand.grid(x = seq(from = Lon.range[1], to = Lon.range[2], by = 0.0001),
                    y = seq(from = Lat.range[1], to = Lat.range[2], by = 0.0001))

coordinates(grid) <- ~x + y #assign coordinates to grid
gridded(grid) <- TRUE ## Create SpatialPixel object

plot(grid, cex = 2, col = "grey")
points(d, pch = 1, col = "blue", cex = 1)
#idw formulae

idw <- idw(formula = sal0 ~ 1, locations = d, newdata = grid)  
idwO = as.data.frame(idw)  
names(idwO)
names(idwO)[1:3] <- c("x", "y", "var1.pred")

ggplot() + geom_tile(data = idwO, aes(x = x, y = y, fill = var1.pred)) 

ggplot() + geom_tile(data = idwO, aes(x = x, y = y, fill = var1.pred))+ scale_fill_gradient(low = "blue", high = "orange")

viet = shapefile("C:/Users/eliecer.diaz/Documents/Recruitment/iho/FinalMap.shp")
crs(viet) = CRS('+init=EPSG:4326')
plot(viet, axes=T)

vietC <- fortify(viet)
##fortify before
##displaying shapefile in ggplot2, i.i
#convert map data to data frame

ggplot() + geom_tile(data = idwO, aes(x = x, y = y, fill = var1.pred)) +
    scale_fill_gradient(low = "blue", high = "orange") + 
    geom_path(data = vietC, aes(long, lat, group = group), colour = "black")


#2)IDW Using the Vornoi method-provide user defined grid
d = DA04
library(dismo)
coordinates(d) <- ~ Lon + Lat
crs(v) <- CRS('+init=EPSG:4326')
v <- voronoi(d)

plot(v)
#summarizes spatial variables
va <- aggregate(viet) #sp package
#set boundaries to vietnam
vca <- intersect(v, va)
spplot(vca, 'sal0', col.regions=rev(get_col_regions()))

#build a raster of vietnam with stated resolution
r <- raster(va, res=0.0001) #1 degree=111 sq km
crs(r) <- CRS('+init=EPSG:4326')
plot(vca)
#projection(r)=CRS("+proj=longlat +ellps=WGS84")

#rasterize polygon
vr <- rasterize(vca, r, 'sal0')
plot(vr)

###############################################
# contour 4 ---Kriging (interpolation)
###############################################

library(gstat)
d = DA04
coordinates(d) <- ~Lon + Lat
summary(d)

Lon.range <- as.numeric( c(23.1372, 23.3031))
Lat.range <- as.numeric( c(59.7575, 59.8640))

grid <- expand.grid(x = seq(from = Lon.range[1], to = Lon.range[2], by = 0.0001),
                    y = seq(from = Lat.range[1], to = Lat.range[2], by = 0.0001))

coordinates(grid) <- ~x + y #assign coordinates to grid
gridded(grid) <- TRUE ## Create SpatialPixel object
crs(d) <- CRS('+init=EPSG:4326')

plot(grid, cex = 2, col = "grey")
points(d, pch = 1, col = "blue", cex = 1)
names(d)
crs(grid) <- CRS('+init=EPSG:4326')
d <- d[!duplicated(d$site),]
v = variogram(log(temp0)~1, d) 
plot(v) #plot semi-variogram
v.fit = fit.variogram(v, vgm("Log")) 
krigeM <- krige(log(temp0) ~ 1, d, grid, model=v.fit)

#display
krigeM %>% as.data.frame %>% ggplot(aes(x=x, y=y)) +
         geom_tile(aes(fill=var1.pred)) + coord_equal() +
                scale_fill_gradient(low = "yellow", high="red") +theme_bw()














##########################################################################
# Map shit!
glgmap   <- get_map(location = c(23.107873, 59.769191, 23.290117, 59.874097),
                    zoom = 11, #17,
                    maptype= "terrain",
                    col = "bw") 

p <- ggmap(glgmap)

p <- p + geom_point(aes(x = Lon,y= Lat),
                        color = "red",
                      data = MA1[MA1$year == 1,]) 

p <- p + geom_text(aes(x = Lon, y = Lat, colour="black", label = site),
                  data = MA1[MA1$year == 1,],
                  size=2,
                  color = "black") 
p


# Interactive
library(plotly)
p <- plot_ly(DA1, x = ~Lat, y = ~Lon, z = ~zerofive, color = ~year,
             marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 150)) 
p %>% layout(title = 'Recruitment variation per year',
         scene = list(xaxis = list(title = 'Longitude',
                      gridcolor = 'rgb(255, 255, 255)',
                      #range = c(2.003297660701705, 5.191505530708712),
                      #type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
               yaxis = list(title = 'Latitude',
                      gridcolor = 'rgb(255, 255, 255)',
                      #range = c(36.12621671352166, 91.72921793264332),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
               zaxis = list(title = 'Recruits',
                            gridcolor = 'rgb(255, 255, 255)',
                            #type = 'log',
                            zerolinewidth = 1,
                            ticklen = 5,
                            gridwith = 2)),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="scatter3d-bubble")
chart_link
#########################################
#############################################
#Data exploration
# A Outliers in Y / Outliers in X
# ----NO
# B Collinearity X
# ----Yes
# C Relationships Y vs X
# ----NO
# D Spatial/temporal aspects of sampling design (not relevant here)
# E Interactions (is the quality of the data good enough to 
#                 include them?)
#F Zero inflation Y
#G Are categorical covariates balanced?
######################################################
source("C:/Users/eliecer.diaz/Documents/RCourseInla/AllRCode/HighstatLibV10.R")
corvif(MA2[,MyVar])
names(MA2)
str(MA2)
#########################################################
###########Fitting the model
#########################################################
M0 <- lm(zerofive ~ Isaeus + sal0 + temp0 + site +
              year, data = MA2)
summary(M0)

M01 <- lmer(zerofive ~   Isaeus + sal0 + temp0 + year +
             (1| site), 
           data = MA2)

AIC(M, M01)

summary(M001)
M <- M01
E0 <- resid(M) #, type = "pearson")
F0 <- fitted(M)
plot(F0, E0)
abline(h = 0, lty = 2)

N  <- nrow(MA2)
p  <- length(fixef(M)) + 2   # 3 sigmas
sum(E0^2) / (N - p)
#0.0257746

MA2$zerofive
library(gamm)
GAM <- gamm(zerofive ~  temp0,# bs = "cr"),
                         data = MA2,
                         family = neg.bin(theta = c(23.64)), #Gamma(link = log),
                         random = list(site =~ 1))

M4 <- glmmadmb(H ~  Region  + (1| station),  
               data = MA2,
               family= "Gamma")



AIC(M1,M2, M3,M3_1$lme, M4) # Exp out!
df      AIC
M1  8 179.5204
M2 10 189.0880
M3  6 176.5584
M4  5 174.5912
summary(M4)
M1 <- M3
E1 <- resid(M1, type = "pearson")
N  <- nrow(MA2)
p  <- length(fixef(M1)) + 2    # 3 sigmas
sum(E1^2) / (N - p)

F1 <- fitted(M1)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values",
     ylab = "Pearson residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)


# Patterns ###NO se que significa!
par(mfrow = c(1,1))
plot(MA2$H, as.numeric(F1), col = MA2$Region)
abline(coef = c(0, 1), lty = 2)   
#########################no good fit, GLM POisson then is not good idea!!

####to test there is non-linearity use GAM

T1 = gam(E1~ s(F1), data =MA2)
plot(T1)

summary(T1) # non significant so non-linearity is confirmed
T2 = gam(E1~s(Bmass), data =MA2)
summary(T2) # non. significant neither
plot(T2)
summary(T2)


levels(MA2$Region)[1] <- c("Hang�")
levels(MA2$Region)[2] <- c("Tv�rminne")
levels(MA2$Region)[3] <- c("Jussar�")

#############################################
##################################PLOTTING
#############################################
library(plyr)
library(ggplot2)

range(MA2$H)
range(MA2$Bmass)
summary(M1)
summary(M3_1$gam)
head(MyData)
MyData <- ddply(MA2, .(Region), summarize,
                Bmass = seq(min(Bmass), 5, length = 30))

P <- predict(M3, newdata = MyData, se = TRUE) #, se = TRUE, type = "response"

MyData$mu <- exp(P$fit )
MyData$ub <- exp(P$fit + 1.96 * P$se.fit )
MyData$lb <- exp(P$fit - 1.96 * P$se.fit )

p <- ggplot()
'''
p <- p + geom_point(data = MA2, 
                    aes(y =H, x = Bmass),
                    shape = 16, 
                    size = 3)
'''
p <- p + xlab("Standarized (Mytilus and Fucus biomass)") + ylab("H Shannon")
p <- p + theme(text = element_text(size=15)) + theme_bw()

p <- p + geom_line(data = MyData, 
                   aes(x = Bmass, y = mu), 
                   colour = "black")


p <- p + geom_ribbon(data = MyData, 
                     aes(x = Bmass, 
                         ymax = ub, 
                         ymin = lb ),
                     alpha = 0.5)

p <- p + facet_grid( ~ Region, scales = "fixed")
p

plot(GAM$gam, select = 1)


    © 2019 GitHub, Inc.
    Terms
    Privacy
    Security
    Status
    Help

    Contact GitHub
    Pricing
    API
    Training
    Blog
    About

