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



known.model$par
neglogpseudolik.config(param = NULL,
                       config = X,
                       crf = known.model,
                       cond.en.form = "feature",
                       ff = f0)

symbolic.conditional.energy(config = c(2,2,1,1,2), condition.element.number = 3, crf = known.model, ff = f0, printQ = F)
symbolic.conditional.energy(config = complement.at.idx(configuration = c(2,2,1,1,2), complement.index = 3), condition.element.number = 3, crf = known.model, ff = f0, printQ = F)



