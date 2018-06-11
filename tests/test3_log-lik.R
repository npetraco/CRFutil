# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)


# Graph formula:
grphf <- ~A:B + A:C + B:C

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj
n.states <- 2

mrf.fit <- make.crf(adj, n.states)
mrf.fit <- make.features(mrf.fit)
mrf.fit <- make.par(mrf.fit, 6)

mrf.fit$node.par[1,1,1]    <- 1
mrf.fit$node.par[2,1,1]    <- 2
mrf.fit$node.par[3,1,1]    <- 3

mrf.fit$edges
mrf.fit$edge.par[[1]][1,1,1] <- 4
mrf.fit$edge.par[[1]][2,2,1] <- 4
mrf.fit$edge.par[[2]][1,1,1] <- 5
mrf.fit$edge.par[[2]][2,2,1] <- 5
mrf.fit$edge.par[[3]][1,1,1] <- 6
mrf.fit$edge.par[[3]][2,2,1] <- 6

mrf.fit$node.par
mrf.fit$edge.par
mrf.fit$n.states
mrf.fit$adj.nodes
mrf.fit$edges


fake.sample <- rbind(
  c(1,1,1),
  c(1,2,1),
  c(1,1,1)
)

UGM_MRF_computeSuffStat(mrf.fit, instances = fake.sample)
mrf.stat(mrf.fit, fake.sample)
