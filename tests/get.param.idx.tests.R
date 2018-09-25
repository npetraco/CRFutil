library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")

n.states <- 2
fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 2)
fitc$node.par[1,1,] <- 1
fitc$node.par[2,1,] <- 1
fitc$edge.par[[1]][1,1,1] <- 2
fitc$edge.par[[1]][2,2,1] <- 2

# Node numbers
nod.idxs <- c(1,2)
# Possible state configurations:
configs <- rbind(
  c(1,1),
  c(1,2),
  c(2,1),
  c(2,2)
)

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

phi.component(config = c(1,2), i = 1, node.par = fitc$node.par, ff = f0)
phi.component(config = c(1,2), i = 1, j=2,
              #node.par = fitc$node.par,
              edge.par = fitc$edge.par,
              edge.mat = fitc$edges, ff = f0)


get.par.idx(config = c(1,2), i = 1, j=2, edge.par = fitc$edge.par, edge.mat = fitc$edges, ff = f0)
get.par.idx(config = c(1,2), i = 1, node.par = fitc$node.par, edge.mat = fitc$edges, ff = f0)
get.node.idxs(par.idx = 1, node.par = fitc$node.par, edge.mat = fitc$edges)
get.node.idxs(par.idx = 2, node.par = fitc$node.par, edge.par = fitc$edge.par, edge.mat = fitc$edges)


