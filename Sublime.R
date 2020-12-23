MyVar <- c("adul_juv", "Biomass", "sal0", "temp0", "four", "two")

Mydotplot(MA1[,MyVar])

#Pairs
Mypairs(MA1[,MyVar])
corvif(MA1[,MyVar]

MA1 <- MA1[!MA1$two > 1500,]
