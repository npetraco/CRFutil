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

# Make up a test parameter vector:
theta <- c(1,2,3,4,5,6)
#theta <- runif(6)
fitc$par <- theta
jk <- make.pots(parms=fitc$par, crf = fitc, rescaleQ=F, replaceQ=T, format = "regular", printQ=T)

fitc$node.pot
fitc$edge.pot

fitc$node.par
fitc$edge.par

# Check: get the parameter vector back?
theta
make.par.from.potentials(fitc)

#
for(xn in 1:nrow(configs)) {
  for(cn in 1:3){
    cond.en <- conditional.config.energy2(config = configs[xn,], condition.element.number = cn, crf=fitc, ff=f0, printQ=FALSE)
    print(paste0("X=(",configs[xn,1], configs[xn,2], configs[xn,3], ") ", "i=", cn, " E(X_i|X/X_i)=",cond.en))
  }
}



conditional.config.energy( config = configs[xn,], condition.element.number = cn, crf=fitc, ff=f0, printQ=FALSE)


