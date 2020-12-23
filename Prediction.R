fGam22 = y ~ -1 + Intercept + fYear2005 + fYear2006 + fYear2007 + Temp1 +
    Temp2 + Temp3 + Temp4 + Sal1 + Sal2 + Sal3 + Sal4 + Adul_juv1 +
    Adul_juv2 + Adul_juv3 + Adul_juv4 + Biom1 + Biom2 + Biom3 +
    Biom4 + f(w, model = spdePois)

range(MA1$sal0)     # 5.203619 6.201219     
range(MA1$temp0)    # 14.93051 18.69009
range(MA1$Biomass)  # 0.000000 8.705393
range(MA1$adul_juv) # 0  2805

range(MA1$sal0.std)     # -1.697723  1.806573     
range(MA1$temp0.std)    # -1.852184  1.929908
range(MA1$Biom.std)     # -1.852184  1.929908
range(MA1$adul_juv.std) # -1.304241  3.426245

MyPData <- expand.grid( 
           Biomass    = seq(1.47, 2.94, length = 25),
           adul_juv       = seq(6.3, 74.25, length = 25),
           fForested = levels(iph$fForested))

Covariates <- data.frame(
  Intercept   = rep(1, N),
  fYear2005   = XYear[,"year2005"],
  fYear2006   = XYear[,"year2006"],
  fYear2007   = XYear[,"year2007"],
  Biom.std = MA1$Biom.std,
  Sal1 = BasisSal[,"Sal1"],
  Sal2 = BasisSal[,"Sal2"],
  Sal3 = BasisSal[,"Sal3"],
  Sal4 = BasisSal[,"Sal4"],
  Temp1 = BasisTemp[,"Temp1"],
  Temp2 = BasisTemp[,"Temp2"],
  Temp3 = BasisTemp[,"Temp3"],
  Temp4 = BasisTemp[,"Temp4"],
  Adul_juv1 = BasisAdul_juv[,"Adul_juv1"],
  Adul_juv2 = BasisAdul_juv[,"Adul_juv2"],
  Adul_juv3 = BasisAdul_juv[,"Adul_juv3"],
  Adul_juv4 = BasisAdul_juv[,"Adul_juv4"],
  Biom1 = BasisBiom[,"Biom1"],
  Biom2 = BasisBiom[,"Biom2"],
  Biom3 = BasisBiom[,"Biom3"],
  Biom4 = BasisBiom[,"Biom4"]
    
  )

  dim(Covariates)
  View(Covariates)

  extractSmoother(G2NB)

##############################################
  ############# PREDICTION #############
##############################################
library(mgcv)

BasisSal <- smoothCon(s(sal0.std, k = 5, fx = TRUE), 
                       data = MA1, 
                       knots = NULL, absorb.cons = TRUE)[[1]]$X

BasisTemp <- smoothCon(s(temp0.std, k = 5, fx = TRUE), 
                       data = MA1, 
                       knots = NULL, absorb.cons = TRUE)[[1]]$X

BasisAdul_juv <- smoothCon(s(adul_juv.std, k = 5, fx = TRUE), 
                       data = MA1, 
                       knots = NULL, absorb.cons = TRUE)[[1]]$X

BasisBiom <- smoothCon(s(Biom.std, k = 5, fx = TRUE), 
                       data = MA1, 
                       knots = NULL, absorb.cons = TRUE)[[1]]$X

colnames(BasisSal) <- paste("Sal", 1:4, sep = "")
colnames(BasisTemp) <- paste("Temp", 1:4, sep = "")
colnames(BasisAdul_juv) <- paste("Adul_juv", 1:4, sep = "")
colnames(BasisBiom) <- paste("Biom", 1:4, sep = "")

library(plyr)
le = 25
Covariates_P <- ddply(Covariates, .(year), summarize,
    Sal1 = seq(-1.063523,  1.223999, length = le),
    Sal2 = seq(-0.6827941,  0.7347652, length = le),
    Sal3 = seq(-0.7904723,  0.4708276, length = le),
    Sal4 = seq(-1.699156,  1.808099, length = le),
    Temp1 = seq(-1.076855,  1.329917, length = le),
  Temp2 = seq(-0.8249484,  0.7065871, length = le),
  Temp3 = seq(-0.8904282,  0.4410686, length = le),
  Temp4 = seq(-1.853748,  1.931538, length = le),
  Adul_juv1 = seq(-0.7780244,  3.6115340, length = le),
  Adul_juv2 = seq(-1.2883198,  0.8357809, length = le),
  Adul_juv3 = seq(-1.0083454,  0.4541147, length = le),
  Adul_juv4 = seq(-1.305342,  3.429137, length = le),
  Biom1 = seq(-0.7461399,  3.8244655, length = le),
  Biom2 = seq(-1.0506059,  0.9236841, length = le),
  Biom3 = seq(-1.2076711,  0.4439852, length = le),
  Biom4 = seq(-1.248955,  3.557343, length = le)

)
Covariates_P$Intercept = 1                
zerofive = seq(0, 2976, length = le)

#Basis smoother NO INTERACTION YET
# lcs.Sal <- inla.make.lincombs(Sal1 = BasisSal[,"Sal1"],
#                                 Sal2 = BasisSal[,"Sal2"],
#                                 Sal3 = BasisSal[,"Sal3"],
#                                 Sal4 = BasisSal[,"Sal4"])

# lcs.Temp <- inla.make.lincombs(Temp1 = BasisTemp[,"Temp1"],
#                                 Temp2 = BasisTemp[,"Temp2"],
#                                 Temp3 = BasisTemp[,"Temp3"],
#                                 Temp4 = BasisTemp[,"Temp4"])

# lcs.Adul_juv <- inla.make.lincombs(Adul_juv1 = BasisAdul_juv[,"Adul_juv1"],
#                                     Adul_juv2 = BasisAdul_juv[,"Adul_juv2"],
#                                     Adul_juv3 = BasisAdul_juv[,"Adul_juv3"],
#                                     Adul_juv4 = BasisAdul_juv[,"Adul_juv4"])

# lcs.Biom <- inla.make.lincombs(Biom1 = BasisBiom[,"Biom1"],
#                                     Biom2 = BasisBiom[,"Biom2"],
#                                     Biom3 = BasisBiom[,"Biom3"],
#                                     Biom4 = BasisBiom[,"Biom4"])


# names(lcs.Sal)   <- paste(names(lcs.Sal), "Sal", sep = "")
# names(lcs.Temp)   <- paste(names(lcs.Temp), "Temp", sep = "")
# names(lcs.Adul_juv)   <- paste(names(lcs.Adul_juv), "Adul_juv", sep = "")
# names(lcs.Biom)   <- paste(names(lcs.Biom), "Biom", sep = "")

# All.lcs <- c(lcs.Sal, lcs.Temp, lcs.Adul_juv, lcs.Biom) # PUTTING TOGETHER

N <- nrow(Covariates_P)

Stack_Pred <- inla.stack(
  tag = "Predict",
  data = list(y = NA),  
  A = list(1, 1),                 
  effects = list(
    Intercept = rep(1, nrow(Covariates_P)),
    Covariates = Covariates_P))

# Combine stacks
All.stacks <- inla.stack(Stack, Stack_Pred)

fGam22 = y ~ -1 + Intercept + fYear2005 + fYear2006 + fYear2007 + Temp1 +
    Temp2 + Temp3 + Temp4 + Sal1 + Sal2 + Sal3 + Sal4 + Adul_juv1 +
    Adul_juv2 + Adul_juv3 + Adul_juv4 + Biom1 + Biom2 + Biom3 +
    Biom4 + f(w, model = spdePois)

G2NB_P <- inla(fGam22,
           family = "nbinomial", 
           #lincomb = All.lcs, #<---HERE THERE ARE THE MF Smoothers
           #control.inla = list(strategy = "gaussian"),
           data = inla.stack.data(All.stacks),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(All.stacks)))

index.Pred <- inla.stack.index(All.stacks,
                              tag = "Predict")$data
F <- G2NB$summary.fitted.values[index.Pred, c("mean", "0.025quant", "0.975quant")]  #210 by 3

# extract smoother
# extract fitted values
# plot shit!
X <- as.matrix(Covariates_P) # ok
beta <- G2NB_P$summary.fixed[,"mean"] # yes
wpm <- G2NB_P$summary.random$w$mean # yes
eta <- X %*% beta #+ A.PoisRepl %*% wpm # yes
mu <- exp(eta) # 
ExpY <- mu / (1 - exp(-mu))
VarY <- ExpY * (1 + mu - ExpY)
Ezt <- (SA$A_tobianus - ExpY) / sqrt(VarY)
Ezt <- as.vector(Ezt)
SA$Ezt <- as.vector(Ezt)
