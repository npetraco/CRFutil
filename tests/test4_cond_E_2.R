# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)

# Check and see if Pr(X) ~= prod Pr(Xi|X/Xi)

# Graph formula: Schmidt Chain
grphf <- ~A:B + B:C + C:D

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj
n.states <- 2

mc <- make.crf(adj, n.states)
# These are what CRF takes as inputs/fits. The "potentials".
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi23 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi34 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))


mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi23
mc$edge.pot[[3]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot

s1<-1
s2<-2
st.sp <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
st.sp
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# Check and see if Pr(X) ~= prod Pr(Xi|X/Xi)
en.dist.info <- distribution.from.energies(st.sp, mc$edges, pot.info$node.energies, pot.info$edge.energies, config.energy, f0)
en.dist <- cbind(st.sp,en.dist.info$state.probs)
en.dist

st.idx <- 5
st.sp[st.idx,]
en.dist[st.idx,]

# Conditional energies??
ce <- sapply(1:4,function(xx){conditional.config.energy(st.sp[st.idx,],
                                                  condition.element.number=xx,
                                                  adj.node.list = mc$adj.nodes,
                                                  edge.mat = mc$edges,
                                                  one.lgp = pot.info$node.energies,
                                                  two.lgp = pot.info$edge.energies,
                                                  ff = f0)})

# Complement conditional energoes??
cce <- sapply(1:4,function(xx){conditional.config.energy(complement.at.idx(st.sp[st.idx,],xx),
                                                  condition.element.number=xx,
                                                  adj.node.list = mc$adj.nodes,
                                                  edge.mat = mc$edges,
                                                  one.lgp = pot.info$node.energies,
                                                  two.lgp = pot.info$edge.energies,
                                                  ff = f0)})

exp(ce)
exp(cce)


pr.ce  <- exp(ce)/(exp(ce) + exp(cce))
pr.cce <- exp(cce)/(exp(ce) + exp(cce))
pr.ce
pr.cce
pr.ce + pr.cce

en.dist[st.idx,]
prod(pr.ce)

