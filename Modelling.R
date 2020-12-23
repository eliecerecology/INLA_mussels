library(haven)
library(tidyverse)
# https://willhipson.netlify.app/post/bayesian_mlm/bayesian_mlm/
#LOO
#https://discourse.mc-stan.org/t/using-loo-with-brms/6736/2
#PREDICTING
#https://rdrr.io/cran/brms/man/predict.brmsfit.html

curran_dat <- read_sav("/home/elvi/CurranLong.sav") %>%
  dplyr::select(id, occasion, read, homecog) %>%
  filter(complete.cases(.))
str(data.frame(curran_dat))

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

library('maps')
library('maptools')
library('mapdata')
data("worldHiresMapEnv")
source("/home/elvi/Documents/MegaSync/RCourseInla/AllRCode/HighstatLibV11.R")
source("spde-tutorial-functions.R")
load(file = "/home/elvi/Documents/MegaSync/Recruitment/MA1.Rda")
names(MA1)
MA1$site <- factor(as.numeric(MA1$site))

LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

xy <- LongLatToUTM(x = MA1$Lon, y = MA1$Lat, zone = 35)
MA1$Xutm <- xy[,2]
MA1$Yutm <- xy[,3]

MA1$Xkm <- MA1$Xutm / 1000
MA1$Ykm <- MA1$Yutm / 1000

# Searching for outliers
MyVar <- c("nine", "four","two", "one", "zerofive", 
             "Biomass", "sal0", "temp0",
             "Xkm", "Ykm", "adul_juv")


Mydotplot(MA1[,MyVar])


#Pairs
Mypairs(MA1[,MyVar])
library(ggcorrplot)
corr <- round(cor(MA1[,MyVar]), 1)
ggcorrplot(corr, hc.order = TRUE, type = "upper",
   lab = TRUE)

# REMOVE BIOMASS
# REMOVE TEMP0
MyVar <- c("nine", "four","two", "one", "zerofive", 
             "sal0", "Xkm", "Ykm")

corvif(MA1[,MyVar])

#MA1$Isaeus.std <- MyStd(MA1$Isaeus)
MA1$sal0.std <- MyStd(MA1$sal0)
#MA1$temp0.std <- MyStd(MA1$temp0)
#MA1$Biom.std <- MyStd(MA1$Biomass)
MA1$nine.std <- MyStd(MA1$nine)
MA1$four.std <- MyStd(MA1$four)
MA1$two.std <- MyStd(MA1$two)
MA1$one.std <- MyStd(MA1$one)
MA1$zerofive.std <- MyStd(MA1$zerofive)
#MA1$adul_juv.std <- MyStd(MA1$adul_juv)


MyVarX <- c("zerofive.std",
            "nine.std",
            "four.std",
            "two.std",
            "one.std",
            "Xkm", 
            "Ykm")

corvif(MA1[,MyVarX])

MyMultipanel.ggp2(MA1,
                    varx = MyVarX,
                    vary = "recr_aut",
                    ylab = "recr_aut",
                    addSmoother = FALSE,
                    addRegressionLine = TRUE,
                    addHorizontalLine = TRUE)

MA2 <- MA1[,MyVarX]
dim(MA2)
names(MA2)
names(MA1)
#############################
##### ANALYSES
#############################

library(brms) # for the analysis
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)

model1 <- brm(bf(zerofive ~  nine.std + four.std + two.std + one.std + sal0.std + (1|site)),
                  family = gaussian,
                  prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 1), class = b),
                       prior(cauchy(0, 1), class = sd),
                       prior(cauchy(0, 1), class = sigma)),
                              data = MA1,
                              iter = 2000,
                              chain = 4,
                              cores = 3)

summary(model1)
summary(model1)$fixed
summary(model1)$random
plot(model1)

set.seed(25)

ids <- sample(unique(MA1$site), 30)
zerofive <- MA1$zerofive
MA1$id <- MA1$site

MA1 %>%
  bind_cols(as_tibble(fitted(model1))) %>%
  filter(id %in% ids) %>%
  ggplot() +
  geom_point(aes(x = sal0.std, y = zerofive), size = 4, alpha = .75, color = "dodgerblue2") +
  geom_point(aes(x = sal0.std, y = Estimate), shape = 1, size = 4, stroke = 1.5) +
  labs(x = "Sites",
       y = "baby mussels",
       title = "Model 1: One Random Effect, 5 Covariates",
       subtitle = "Blue points are observed values. Black circles are fitted values.") +
  scale_x_continuous(expand = c(.075, .075), breaks = 0:3) +
  facet_wrap(~id, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = .5))

print(model1)
class(MA1$site)
class(MA1$sal0.std)
MA1$sal0.std <- factor(MA1$sal0.std)
##########
# Model 2
##########

model2 <- brm(bf(zerofive ~  1 + (1|sal0.std) + (1|site)),
                  family = gaussian,
                  prior = c(prior(normal(0, 10), class = Intercept),
                       #prior(normal(0, 1), class = b),
                       prior(cauchy(0, 1), class = sd),
                       prior(cauchy(0, 1), class = sigma)),
                              data = MA1,
                              iter = 4000,
                              chain = 4,
                              cores = 3)
summary(model2)
MA1 %>%
  bind_cols(as_tibble(fitted(model2))) %>%
  filter(id %in% ids) %>%
  ggplot() +
  geom_point(aes(x = four.std, y = zerofive), size = 4, alpha = .75, color = "dodgerblue2") +
  geom_point(aes(x = four.std, y = Estimate), shape = 1, size = 4, stroke = 1.5) +
  labs(x = "Sites",
       y = "baby mussels",
       title = "Model 1: One Random Effect, 5 Covariates",
       subtitle = "Blue points are observed values. Black circles are fitted values.") +
  scale_x_continuous(expand = c(.075, .075), breaks = 0:3) +
  facet_wrap(~id, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = .5))

model3 <- brm(bf(zerofive ~  1 + (1|sal0.std) + (1|site)),
                  family = poisson,
                  prior = c(prior(normal(0, 10), class = Intercept)),
                       #prior(normal(0, 1), class = b),
                       #prior(cauchy(0, 1), class = sd),
                       #prior(normal(0, 1), class = sigma)),
                              data = MA1,
                              iter = 4000,
                              chain = 4,
                              cores = 3)
summary(model2)
MA1 %>%
  bind_cols(as_tibble(fitted(model3))) %>%
  filter(id %in% ids) %>%
  ggplot() +
  geom_point(aes(x = four.std, y = zerofive), size = 4, alpha = .75, color = "dodgerblue2") +
  geom_point(aes(x = four.std, y = Estimate), shape = 1, size = 4, stroke = 1.5) +
  labs(x = "Sites",
       y = "baby mussels",
       title = "Model 1: One Random Effect, 5 Covariates",
       subtitle = "Blue points are observed values. Black circles are fitted values.") +
  scale_x_continuous(expand = c(.075, .075), breaks = 0:3) +
  facet_wrap(~id, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = .5))

ggplot(data = MA1, 
       aes(x   =four.std,
           y   = zerofive,
           col = as.factor(site)))+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_point(size     = .7,
             alpha    = .8,
             position = "jitter")+
  geom_smooth(method = lm,
              se     = FALSE, 
              size   = 2,
              alpha  = .8)+
  theme_minimal()+
  labs(title    = "Linear Relationship for Different Years of Teacher Experience as Observed", 
       subtitle = "The linear relationship between the two is not the same for all classes", 
       col      = "Years of\nTeacher\nExperience")

MA1 %>%
  ggplot(aes(x = sal0.std, y = zerofive, group = site)) +
  geom_line(size = .75, alpha = .20) +
  labs(x = "Assessment Period",
       y = "Reading Ability") +
  theme_minimal(base_size = 16)


#MODEL 4
model4 <- brm(bf(recr_aut ~ 1 + zerofive.std + (1|year) + (1|site)),
                  family = poisson,
                  prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 1), class = b)),
                       #prior(cauchy(0, 1), class = sd),
                       #prior(normal(0, 1), class = sigma)),
                              data = MA1,
                              iter = 4000,
                              chain = 4,
                              cores = 4)

summary(model4)

iys <- sample(unique(MA1$year), 4)
zerofive <- MA1$zerofive
MA1$iy <- MA1$year


MA1 %>%
  bind_cols(as_tibble(fitted(model4))) %>%
  filter(iy %in% iys) %>%
  ggplot() +
  geom_point(aes(x = zerofive.std, y = recr_aut), size = 4, alpha = .75, color = "dodgerblue2") +
  geom_point(aes(x = zerofive.std, y = Estimate), shape = 1, size = 4, stroke = 1.5) +
  labs(x = "Sites",
       y = "mussels in ropes",
       title = "Model 1: One Random Effect, 1 Covariates",
       subtitle = "Blue points are observed values. Black circles are fitted values.") +
  scale_x_continuous(expand = c(.075, .075), breaks = 0:3) +
  facet_wrap(~iy, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = .5))


fit_nb <- brm(bf(Effort  ~ 1 + Input + Output + Enquiry + File + 
                   Interface + Added + Changed + (1 | DQR)),
              family = negbinomial,
              prior = c(prior(normal(0,1), class=b),
                        prior(normal(0,1), class=Intercept)),
              cores=4, chains=4, data = dataset, 
              control = list(adapt_delta=0.99, max_treedepth=13)
)
###################
###### STAN  ######
###################

summary(model1)
model1tranformed <- ggs(model1) # 

library(rstan)  
mod <- "
data {
  int<lower=1> N; //the number of observations
  int<lower=1> K; //number of predictiors
  matrix[N,K] x; //predictor matrix
  vector[N] y; //the response variable
}
parameters {
  real alpha; //intercept         
  vector[K] beta; //matrix of group-level regression coefficients
  real<lower=0> sigma; //standard deviation of the individual observations
}
model {
  y ~ normal(alpha + x*beta, sigma); // likelihood
}

"
View(MA2)

N = nrow(MA2)
x <- MA2
K <- ncol(MA2)
y = MA1$zerofive
data_stan <- list(y = y,
                 x = x,
                 N= N,
                 K =K)

rstan_options(auto_write = FALSE)
options(mc.cores = 3)
getwd()
mod2 <- stan(model_code = mod,
                            data=data_stan,
                            iter = 10000, warmup = 2000)

traceplot(mod2)
