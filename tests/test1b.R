library(CRFutil)
library(Rgraphviz)

# Graph formula:
grphf <- ~A:B + B:C + C:D + D:A

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj
n.states <- 2

samp.info  <- sim.field.random(adj, n.states, 100)
true.model <- samp.info$model
true.model$node.pot
true.model$edge.pot


ft.ifo0 <- mrf.standard.fit(samp.info$samples, grphf, n.states, mrf.exact.nll, infer.exact)
ft.ifo1 <- mrf.standard.fit(samp.info$samples, grphf, n.states, mrf.junction.nll, infer.junction)

ft.ifo0$inference.info$logZ
ft.ifo1$inference.info$logZ

pot.info0 <- make.gRbase.potentials(ft.ifo0$fit.model, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info0

pot.info1 <- make.gRbase.potentials(ft.ifo1$fit.model, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info1

ft.ifo0$fit.model$edges
ft.ifo1$fit.model$edges

bel.info0 <- make.gRbase.beliefs(ft.ifo0$inference.info, node.names = gp@nodes, edge.mat = ft.ifo0$fit.model$edges, state.nmes=NULL)
bel.info1 <- make.gRbase.beliefs(ft.ifo1$inference.info, node.names = gp@nodes, edge.mat = ft.ifo1$fit.model$edges, state.nmes=NULL)
bel.info0$node.beliefs[[1]]
bel.info1$node.beliefs[[1]]

bel.info0$node.beliefs[[2]]
bel.info1$node.beliefs[[2]]

bel.info0$node.beliefs[[3]]
bel.info1$node.beliefs[[3]]

bel.info0$node.beliefs[[4]]
bel.info1$node.beliefs[[4]]

bel.info0$edge.beliefs[[1]]
bel.info1$edge.beliefs[[1]]

bel.info0$edge.beliefs[[2]]
bel.info1$edge.beliefs[[2]]

bel.info0$edge.beliefs[[3]]
bel.info1$edge.beliefs[[3]]

bel.info0$edge.beliefs[[4]]
bel.info1$edge.beliefs[[4]]


# Compute joint from potentials with gRbase
gR.dist.info0 <- distribution.from.potentials(pot.info0$node.potentials, pot.info0$edge.potentials)
gR.dist0 <- as.data.frame(as.table(gR.dist.info0$state.probs))
gR.dist0

# Generate state space
s1<- 1
s2<- 2
f0 <- function(y){ as.numeric(c((y==s1),(y==s2)))}
st.sp <- expand.grid(A=c(s1,s2),B=c(s1,s2),C=c(s1,s2),D=c(s1,s2))
st.sp

# Exact
# Compute joint from energies
en.dist.info0 <- distribution.from.energies(st.sp, ft.ifo0$fit.model$edges, pot.info0$node.energies, pot.info0$edge.energies, energy, f0)
en.dist0 <- cbind(st.sp,en.dist.info0$state.probs)
en.dist0

# Check to make sure dists from both gR and en have the same state indices
colnames(gR.dist0)[-5]
# D,C,A,B is the order in gRdist
gR.distQQ0 <- sapply(1:nrow(st.sp), function(xx){gR.dist.info0$state.probs[st.sp[xx,4], st.sp[xx,3], st.sp[xx,1], st.sp[xx,2]]})
gR.distQQ0

cbind(st.sp, en.dist.info0$state.probs, gR.distQQ0, en.dist.info0$state.probs-gR.distQQ0)

# Junction tree
# Compute joint from potentials with gRbase
gR.dist.info1 <- distribution.from.potentials(pot.info1$node.potentials, pot.info1$edge.potentials)
gR.dist1 <- as.data.frame(as.table(gR.dist.info1$state.probs))
gR.dist1
cbind(gR.dist0[,5], gR.dist1[,5], gR.dist0[,5]-gR.dist1[,5])

# Compute joint from energies
en.dist.info1 <- distribution.from.energies(st.sp, ft.ifo1$fit.model$edges, pot.info1$node.energies, pot.info1$edge.energies, energy, f0)
en.dist1 <- cbind(st.sp,en.dist.info1$state.probs)
en.dist1
cbind(en.dist0[,5],en.dist1[,5], en.dist0[,5]-en.dist1[,5]) # Compare exact and JT.... HUH????


# Check to make sure dists from both gR and en have the same state indices
colnames(gR.dist1)[-5]
# D,C,A,B is the order in gRdist
gR.distQQ1 <- sapply(1:nrow(st.sp), function(xx){gR.dist.info1$state.probs[st.sp[xx,4], st.sp[xx,3], st.sp[xx,1], st.sp[xx,2]]})
gR.distQQ1

cbind(st.sp, en.dist.info1$state.probs, gR.distQQ1, en.dist.info1$state.probs-gR.distQQ1)


library(prodlim)
re.arr.idxs <- sapply(1:nrow(st.sp), function(xx){row.match(x = st.sp[xx,] ,table =  gR.dist0[,c(3,4,2,1)])})
re.arr.idxs

cbind(gR.dist0[re.arr.idxs,5], gR.dist1[re.arr.idxs,5], gR.dist0[re.arr.idxs,5]-gR.dist1[re.arr.idxs,5])
cbind(en.dist0[,5],en.dist1[,5], en.dist0[,5]-en.dist1[,5])

# Pr diffrences Exact vs JT. Left is by gR, right is by energy
cbind(gR.dist0[re.arr.idxs,5]-gR.dist1[re.arr.idxs,5], en.dist0[,5]-en.dist1[,5])
