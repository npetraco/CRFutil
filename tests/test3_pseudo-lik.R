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

mc <- make.crf(adj, n.states)
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
#Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  rbind(
    c(30, 5),
    c(1, 10)
  )
Psi23 <-
  rbind(
    c(100, 1),
    c(1, 100)
  )
Psi31 <-
  rbind(
    c(20, 1),
    c(1, 30)
  )

mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
#mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi31
mc$edge.pot[[3]] <- Psi23
#mc$edge.pot[[4]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot


s1<-1
s2<-2
st.sp <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
st.sp
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# Check and see if Pr(X) ~= prod Pr(Xi|X/Xi)
en.dist.info <- distribution.from.energies(st.sp, mc$edges, pot.info$node.energies, pot.info$edge.energies, config.energy, f0)
en.dist <- cbind(st.sp,en.dist.info$state.probs)
en.dist

pl.dist.info <- pseudolikelihoods.from.energies(
  st.sp, adjacent.nodes = mc$adj.nodes,
  edges.mat = mc$edges,
  node.energies = pot.info$node.energies,
  edge.energies = pot.info$edge.energies,
  conditional.energy.func = conditional.config.energy,
  ff = f0)

cbind(100*round(en.dist[,4],5), 100*round(pl.dist.info$pseudo.likelihoods,5))

pl.dist.info
