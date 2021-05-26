library(CRFutil)
library(Rgraphviz)


# Graph formula:
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


#State space:
s1<-1
s2<-2
#f0 <- function(y){ as.numeric(c((y==s1),(y==s2)))}
#st.sp <- expand.grid(A=c(s1,s2),B=c(s1,s2),C=c(s1,s2),D=c(s1,s2))
st.sp <- as.matrix(expand.grid(A=c(s1,s2),B=c(s1,s2),C=c(s1,s2),D=c(s1,s2)))
st.sp

# Decorate potentials/energies with gRbase conventions:
#pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

# Compute product of potential values for a state config acreoos the nodes:
config.potential(config=c(1,1,1,1), edges.mat=mc$edges, node.pots=pot.info$node.potentials, edge.pots=pot.info$edge.potentials)
config.potential(config=c(2,1,2,2), edges.mat=mc$edges, node.pots=pot.info$node.potentials, edge.pots=pot.info$edge.potentials)

prodPots <- sapply(1:nrow(st.sp), function(xx){config.potential(config=st.sp[xx,], edges.mat=mc$edges, node.pots=pot.info$node.potentials, edge.pots=pot.info$edge.potentials)})
prodPots

