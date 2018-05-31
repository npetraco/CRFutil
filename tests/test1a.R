# Initalize an mrf-object:
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

mc <- make.crf(adj, n.states)
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.9,  0.1) *10

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
Psi34 <-
  rbind(
    c(1, 100),
    c(100, 1)
  )
Psi41 <-
  rbind(
    c(100, 1),
    c(1, 100)
  )

mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi41
mc$edge.pot[[3]] <- Psi23
mc$edge.pot[[4]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot


# Energy func that takes one and two body energies
# Set and assign to MyEnergy??
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==s1),(y==s2)))}

st.sp <- expand.grid(A=c(s1,s2),B=c(s1,s2),C=c(s1,s2),D=c(s1,s2))
st.sp

# Decorate potentials/energies with gRbase conventions:
#pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes, state.nmes = c("+1","-1"))
pot.info

# Compute joint from potentials with gRbase
gR.dist.info <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
as.data.frame(as.table(gR.dist.info$state.probs))

# Compute joint from energies
en.dist.info <- distribution.from.energies(st.sp, mc$edges, pot.info$node.energies, pot.info$edge.energies, energy, f0)
cbind(st.sp,en.dist.info$state.probs)


class(gR.dist.info$state.probs)
class(en.dist.info$state.probs)

st.sp
