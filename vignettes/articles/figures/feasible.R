################################################################################
# Plot possible posterior probabilities for MOI = (2, 1)
################################################################################

library(Pv3Rs)

recru <- project2D(c(1,0,0))
relap <- project2D(c(0,1,0))
reinf <- project2D(c(0,0,1))

# Assume MOI = (2, 1)
# 9 graphs compatible with relapse
# 4 graphs compatible with recrudescence
# 2 graphs compatible with reinfection

relap.recru <- project2D(c(9,4,0)/(9+4))
relap.reinf <- project2D(c(0,2,9)/(9+2))
intersect.pt <- project2D(c(9*2,4*2,9*4)/(9*2+4*2+9*4))

png("feasible.png", width = 4, height = 4, units = "in", res = 300)

par(mar=rep(0,4))
V_labels <- c("Recrudescence", "Relapse", "Reinfection")
plot_simplex(v.labels =  V_labels, v.colours=rep("white",3))
segments(recru['x'], recru['y'], relap.reinf['x'], relap.reinf['y'], lty=2)
segments(reinf['x'], reinf['y'], relap.recru['x'], relap.recru['y'], lty=2)

feasible <- rbind(relap, relap.reinf, intersect.pt, relap.recru)
polygon(x=feasible[,"x"], y=feasible[,"y"], border=NA, col=adjustcolor("green", alpha.f = 0.35))
dev.off()
