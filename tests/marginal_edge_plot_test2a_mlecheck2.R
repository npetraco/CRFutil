library(CRFutil)
library(rstanarm)
library(rethinking)

grf.eq <- ~X + Y + X:Y
XY <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)

XY$node.pot[1,] <- c(0.8,0.2)
XY$node.pot[2,] <- c(1-0.37,0.37)
XY$edge.pot[[1]][1,1] <- XY$edge.pot[[1]][2,2] <- 1 # ******* theta true 11,22 is theta-times more likely than 12 or 21
XY$node.pot
XY$edge.pot

num.samps <- 100000
samps     <- sample.exact(XY, num.samps)
mrf.sample.plot(samps)
XY.emp.prs <- marginal.edge.emp.pr(samps) # Check sample margins. Marginals should be node pots and show independence
XY.emp.prs

# Loglinear Poison MLE
llm <- marginal.edge.loglin(samps, conf.level = 0.95)
llm$glm.poi.rescaled.node.pot
llm$glm.poi.rescaled.edge.pot

# Put into MRF obj
grf.eq2          <- ~A + B + A:B
AB               <- make.empty.field(graph.eq = grf.eq2, parameterization.typ = "standard", plotQ = F)
AB$node.pot      <- llm$glm.poi.rescaled.node.pot
AB$edge.pot[[1]] <- llm$glm.poi.rescaled.edge.pot
AB$node.pot
AB$edge.pot
#AB$par <- llm$glm.theta.est
AB$par <- c(log(AB$node.pot)[,1], log(AB$edge.pot[[1]])[1,1]) # thetas from rescaled pots, just in case
AB$par
exp(AB$par)

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
XY.emp.prs$`Pr(X1,X2)` # Check against empirical values
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
XY.emp.prs$`Pr(X1)`                        # Check against empirical values
sum(AB.bels[["Bel(A)"]])

# Pr(B)
AB.bels$`Bel(B)`
AB.bels$"Bel(B)"
AB.bels[["Bel(B)"]]
AB.bx$node.bel[2,]
ar_marg(AB.bels[["Bel(A,B)"]], marg = "B")
colSums(AB.bx$edge.bel[[1]])
colSums(AB.bels[["Bel(A,B)"]])
XY.emp.prs$`Pr(X2)`                        # Check against empirical values
sum(AB.bels[["Bel(B)"]])

# Pr(A|B) = Pr(A,B)/Pr(B)
ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(B)"]])
# Is Pr(A|B) == Pr(A) ?? Independence test
AB.bels[["Bel(A)"]]
rowSums(ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(B)"]]))
XY.emp.prs$`Pr(X1|X2)`                        # Check against empirical values
XY.emp.prs$`Pr(X1)`                           # Check against empirical values

# Pr(B|A) = Pr(A,B)/Pr(A)
ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(A)"]])
# Is Pr(B|A) == Pr(B) ?? Independence test
AB.bels[["Bel(B)"]]
rowSums(ar_div(AB.bels[["Bel(A,B)"]], AB.bels[["Bel(A)"]]))
XY.emp.prs$`Pr(X2|X1)`                        # Check against empirical values
XY.emp.prs$`Pr(X2)`                           # Check against empirical values

