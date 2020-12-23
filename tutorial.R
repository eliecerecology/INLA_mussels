rm(list=ls(all=TRUE))

library(geoR)
library(ggplot2)
#https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html
View(parana)
ggplot(data.frame(cbind(parana$coords, Rainfall = parana$data)))+
  geom_point(aes(east, north, color = Rainfall), size = 2) +
  coord_fixed(ratio = 1) +
  scale_color_gradient(low = "blue", high = "orange") +
  geom_path(data = data.frame(parana$border), aes(east, north)) +
  theme_bw()

library(INLA)
coo <- parana$coords
plot(coo)
summary(dist(coo))


###
# 1
###
mesh <- inla.mesh.2d(
  loc = coo, offset = c(50, 100), # c(inner, outer)
  cutoff = 1, max.edge = c(30, 60)
)
plot(mesh)
points(coo, col = "red")

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

indexs <- inla.spde.make.index("s", spde$n.spde) # number of vertices
A <- inla.spde.make.A(mesh = mesh, loc = coo)

# dim(as.data.frame(parana))
# dim(A)
# nrow(coo)
# head(A)


##
#PREDICTION
##
#First, we construct a grid called coop with 50×50 locations by using expand.grid()
#and combining vectors x and y which contain coordinates in the range of parana$border.

bb <- bbox(parana$border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 50)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 50)
coop <- as.matrix(expand.grid(x, y))
dim(coop)

ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  parana$border[, 1], parana$border[, 2]
)
plot(parana) #$border[, 1], parana$border[, 2])

coop <- coop[which(ind == 1), ]
plot(coop, asp = 1)

# New projection matrix for prediction
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
dim(Ap)    

stk.e <- inla.stack(
  tag = "est",
  data = list(y = parana$data),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)

stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + f(s, model = spde)

res <- inla(formula,
  data = inla.stack.data(stk.full),
  control.predictor = list(
    compute = TRUE,
    A = inla.stack.A(stk.full)
  )
)

# Results
index <- inla.stack.index(stk.full, tag = "pred")$data
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

dpm <- rbind(
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_mean, variable = "pred_mean"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ll, variable = "pred_ll"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ul, variable = "pred_ul"
  )
)
dpm$variable <- as.factor(dpm$variable)

ggplot(dpm) + geom_tile(aes(east, north, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Rainfall",
    low = "blue", high = "orange"
  ) +
  theme_bw()

newloc <- cbind(c(219, 678, 818), c(20, 20, 160))
Aproj <- inla.spde.make.A(mesh, loc = newloc)
Aproj %*% res$summary.random$s$mean

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh,
  xlim = rang[, 1], ylim = rang[, 2],
  dims = c(300, 300)
)
