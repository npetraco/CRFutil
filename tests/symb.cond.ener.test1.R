library(CRFutil)

# Slayer field model:
grphf <- ~1:2+1:3+1:4+1:5+2:3+2:4+2:5+3:4+3:5
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)


# Make up random potentials and return a CRF-object
slay    <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=1, seed=1)
samps       <- slay$samples
known.model <- slay$model
known.model$node.par
    known.model$edges[9,]
known.model$edge.par[[9]]


# Feature function
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

# Grab the first sampled configuration:
X <- samps[1,]
X

#a.  Gradient of E(X3|X/X3)
dE.X3.Xb <- conditional.energy.gradient(config = X, condition.element.number = 3, crf = known.model, ff = f0, printQ=F)
dE.X3.Xb <- dE.X3.Xb$conditional.grad

# dE/dth7:
dE.X3.Xb[7]

#b. Requisite energies, potentials and gradients
E.X3.Xb  <- conditional.config.energy2(theta.par = NULL, config = X, condition.element.number = 3, crf = known.model, ff = f0)
P.X3.Xb  <- exp(E.X3.Xb)              # potential needed
dE.X3.Xb <- conditional.energy.gradient(config = X, condition.element.number = 3, crf = known.model, ff = f0, printQ=F)
dE.X3.Xb <- dE.X3.Xb$conditional.grad # derivative needed

# Complement energies, potentials and gradients
Xc        <- complement.at.idx(configuration = X, complement.index = 3)
E.X3.Xbc  <- conditional.config.energy2(theta.par = NULL, config = Xc, condition.element.number = 3, crf = known.model, ff = f0)
P.X3.Xbc  <- exp(E.X3.Xbc)              # complement potential needed
dE.X3.Xbc <- conditional.energy.gradient(config = Xc, condition.element.number = 3, crf = known.model, ff = f0, printQ=F)
dE.X3.Xbc <- dE.X3.Xbc$conditional.grad # complement derivative needed

# Assemble \frac{\partial}{\partial \theta_7} Z_{X_3|{\bf X}_1\slash X_3}
dZ.X3.Xb.th7 <- P.X3.Xb * dE.X3.Xb[7] + P.X3.Xbc * dE.X3.Xbc[7]
dZ.X3.Xb.th7

# c. First compute the required Z
logZ.X3.Xb <- logsumexp2(c(E.X3.Xb, E.X3.Xbc))
Z.X3.Xb    <- exp(logZ.X3.Xb)

# Compute the required probabilities
Pr.X3.Xb   <- P.X3.Xb/Z.X3.Xb
Pr.X3.Xbc  <- P.X3.Xbc/Z.X3.Xb

# Assemble \frac{\partial}{\partial \theta_7}\log\Big( Z_{X_3|{\bf X}_1\slash X_3} \Big)
Pr.X3.Xb * dE.X3.Xb[7] + Pr.X3.Xbc * dE.X3.Xbc[7]

# Or we can do this:
1/Z.X3.Xb * dZ.X3.Xb.th7

known.model$par

nlpl.info <- neglogpseudolik.config(
  param        = NULL,
  config       = X,
  crf          = known.model,
  cond.en.form = "feature",
  gradQ        = T,
  ff           = f0)
nlpl.info


gnlpl.info <-grad.neglogpseudolik.config(
  param                      = NULL,
  config                     = X,
  crf                        = known.model,
  cond.en.form               = "feature",
  ff                         = f0)

gnlpl.info$grad.neglogpseudolik

gnlpl.info$alpha
gnlpl.info$alpha[7,3]

gnlpl.info$dZ
gnlpl.info$dZ[7,3]

gnlpl.info$Ealpha
gnlpl.info$Ealpha[7,3]

conditional.energy.gradient(config = X,   condition.element.number = 1, crf = known.model, ff = f0)
symbolic.conditional.energy(config = c(2,2,1,1,2), condition.element.number = 3, crf = known.model, ff = f0, printQ = F)
symbolic.conditional.energy(config = complement.at.idx(configuration = c(2,2,1,1,2), complement.index = 3), condition.element.number = 3, crf = known.model, ff = f0, printQ = F)



