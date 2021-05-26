library(CRFutil)

# Graph formula:
grphf <- ~A:B:C
adj   <- ug(grphf, result="matrix")

n.states <- 2
fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 6)

# Parameterization:
fitc$node.par[1,1,] <- 1
fitc$node.par[2,1,] <- 2
fitc$node.par[3,1,] <- 2
fitc$edge.par[[1]][1,1,1] <- 3
fitc$edge.par[[1]][2,2,1] <- 3
fitc$edge.par[[2]][1,1,1] <- 2
fitc$edge.par[[2]][2,2,1] <- 4
fitc$edge.par[[3]][1,1,1] <- 5
fitc$edge.par[[3]][2,2,1] <- 6

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

configs <- rbind(
  c(1,1,1),
  c(1,2,1),
  c(2,2,1),
  c(2,2,2)
)

# Make up a test parameter vector with theta_1 = 1, theta_2 = 2, etc:
theta <- c(1,2,3,4,5,6)

fitc$par <- theta
jk <- make.pots(parms=fitc$par, crf = fitc,
                rescaleQ=F, replaceQ=T, format = "regular", printQ=T)



# Gradient functionality testing:
cfg <- configs[4,]
cfg
conditional.phi.mat.test(config = cfg, condition.element.number = 1, crf = fitc, ff = f0, printQ=T)
#
