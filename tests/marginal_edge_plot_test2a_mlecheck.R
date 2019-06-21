library(CRFutil)
library(rstanarm)
library(rethinking)


# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 100000

# Samples are independent
samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.37)
)

dev.off()
mrf.sample.plot(samps)

# Loglinear Poison MLE
llm <- marginal.edge.loglin(samps, conf.level = 0.95)
llm$glm.poi.rescaled.node.pot
llm$glm.poi.rescaled.edge.pot

# Put into MRF obj
grf.eq           <- ~A + B + A:B
AB               <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)
AB$node.pot      <- llm$glm.poi.rescaled.node.pot
AB$edge.pot[[1]] <- llm$glm.poi.rescaled.edge.pot
AB$node.pot
AB$edge.pot
#AB$par <- llm$glm.theta.est
AB$par <- c(log(AB$node.pot)[,1], log(AB$edge.pot[[1]])[1,1]) # thetas from rescaled pots, just in case
AB$par

# Compute beliefs
AB.bels <- marginal.edge.bels(AB, node.names = c("A","B"), state.names = c("1","2"))
AB.bx   <- infer.exact(AB)

# Checks:
# Pr(A,B)
AB.bels[["Bel(A,B)"]]
AB.bels[[1]]
AB.bels$`Bel(A,B)`
AB.bels$"Bel(A,B)"
AB.bx$edge.bel[[1]]
sum(AB.bels[[1]])
sum(AB.bx$edge.bel[[1]])

# Pr(A)
AB.bels$`Bel(A)`
AB.bels$"Bel(A)"
AB.bels[["Bel(A)"]]
AB.bx$node.bel[1,]
ar_marg(AB.bels[["Bel(A,B)"]], marg = "A") # Marginalize out A from A,B
rowSums(AB.bx$edge.bel[[1]])               # Also marginalize out A from A,B
rowSums(AB.bels[["Bel(A,B)"]])
sum(AB.bels[["Bel(A)"]])

# Pr(B)
AB.bels$`Bel(B)`
AB.bels$"Bel(B)"
AB.bels[["Bel(B)"]]
AB.bx$node.bel[2,]
ar_marg(AB.bels[["Bel(A,B)"]], marg = "B")
colSums(AB.bx$edge.bel[[1]])
colSums(AB.bels[["Bel(A,B)"]])
sum(AB.bels[["Bel(B)"]])

# Pr(A|B) = Pr(A,B)/Pr(B)
ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(B)"]])
# Is Pr(A|B) == Pr(A) ?? Independence test
AB.bels[["Bel(A)"]]
rowSums(ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(B)"]]))

# Pr(B|A) = Pr(A,B)/Pr(A)
ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(A)"]])
# Is Pr(B|A) == Pr(B) ?? Independence test
AB.bels[["Bel(B)"]]
rowSums(ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(A)"]]))

