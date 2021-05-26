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

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

configs <- rbind(
  c(1,1,1),
  c(1,2,1),
  c(2,2,1),
  c(2,2,2)
)
# Node parameter components:
get.par.idx(config = configs[4,], i = 3, node.par = fitc$node.par, ff = f0)

# Edge parameter components:
get.par.idx(config = configs[2,], i = 2, j = 3, edge.par = fitc$edge.par, edge.mat = fitc$edges, ff = f0)

fitc$par <- log(c(1,2,3,4,5,6))
fitc$par
fitc$adj.edges
fitc$edges

conditional.config.energy2(configs[2,], condition.element.number = 3, crf=fitc, ff=f0, printQ=FALSE)
fitc$adj.edges

get.par.idx(config = c(1,2,1),
            i        = 2,
            j        = 3,
            edge.par = fitc$edge.par,
            edge.mat = fitc$edges,
            ff       = f0)
phi.component(config = c(1,2,1),
                       i        = 2,
                       j        = 3,
                       edge.par = fitc$edge.par,
                       edge.mat = fitc$edges,
                       ff       = f0)

