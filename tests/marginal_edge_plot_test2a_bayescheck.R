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
XY.emp.prs <- marginal.edge.emp.pr(samps) # Check sample margins. Marginals should be node pots and show independence
XY.emp.prs

# Now fit a Bayes loglin model instead
marg.eg.info <- marginal.edge.bayes.loglin(samps)
eg.post.pots <- marg.eg.info$rescaled.posterior.pots

# Look at posterior of exp(omega) = pot.omega. It should comfortably cover 1 since the sample is independent
hist(eg.post.pots[,3])
hist(log(eg.post.pots[,3]))

# Insert a sampled posterior pot into an MRF object ana compute beliefs
i <- 1
grf.eq           <- ~A + B + A:B
AB               <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)
AB$par           <- log(eg.post.pots[i,]) # LOOP HERE
out.mkpot        <- make.pots(parms = AB$par, crf = AB, rescaleQ = F, replaceQ = T)
AB$node.pot
AB$edge.pot
eg.post.pots[1,] # Same as abive

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


