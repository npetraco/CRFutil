library(CRFutil)

# First set up a Slayer field with known potentials with which to generate a sample of X:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E
gp    <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj   <- ug(grphf, result="matrix")
adjN

n.states    <- 2
slay        <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=1000, seed=1)

# A sample of X for which we will fit a model and obtain a theta:
samps       <- slay$samples

# # Fit an MRF to the sample with the intention of obtaining a theta
# # Use the standard parameterization:
# fit <- make.crf(adj, n.states)
# fit <- make.features(fit)
# fit <- make.par(fit, 14)
# 
# # One param per node:
# for(i in 1:nrow(fit$node.par)){
#   fit$node.par[i,1,1] <- i
# }
# fit$node.par
# 
# # Edge parameterization:
# fit$edges # Check edge order first!
# fit$edge.par[[1]][1,1,1] <- 6
# fit$edge.par[[1]][2,2,1] <- 6
# fit$edge.par[[2]][1,1,1] <- 7
# fit$edge.par[[2]][2,2,1] <- 7
# fit$edge.par[[3]][1,1,1] <- 8
# fit$edge.par[[3]][2,2,1] <- 8
# fit$edge.par[[4]][1,1,1] <- 9
# fit$edge.par[[4]][2,2,1] <- 9
# fit$edge.par[[5]][1,1,1] <- 10
# fit$edge.par[[5]][2,2,1] <- 10
# fit$edge.par[[6]][1,1,1] <- 11
# fit$edge.par[[6]][2,2,1] <- 11
# fit$edge.par[[7]][1,1,1] <- 12
# fit$edge.par[[7]][2,2,1] <- 12
# fit$edge.par[[8]][1,1,1] <- 13
# fit$edge.par[[8]][2,2,1] <- 13
# fit$edge.par[[9]][1,1,1] <- 14
# fit$edge.par[[9]][2,2,1] <- 14
#
# fit$edge.par



# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }
1111
# Enumerate all the state configurations
all.configs <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(all.configs) <- gp@nodes
all.configs


# Compare config energies with energies computed as theta \dot phi:
# First compute the “features” (phi) for all each possible configuration
M.all  <- compute.model.matrix(all.configs, fit$edges, fit$node.par, fit$edge.par, f)
