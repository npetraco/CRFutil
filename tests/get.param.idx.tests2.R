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
fitc$edge.par[[3]][2,1,1] <- 4
fitc$edges
fitc$node.par
fitc$edge.par

# Node numbers
X <- c(1,2,1)

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

phi.features(
  config    = X,
  edges.mat = fitc$edges,
  node.par  = fitc$node.par,
  edge.par  = fitc$edge.par,
  ff        = f0
)

get.par.idx(config = X, i = 1, j=2, edge.par = fitc$edge.par, edge.mat = fitc$edges, ff = f0)
get.par.idx(config = X, i = 1, node.par = fitc$node.par, edge.mat = fitc$edges, ff = f0)
get.node.idxs(par.idx = 1, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 2, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 3, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 4, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 5, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 6, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)
