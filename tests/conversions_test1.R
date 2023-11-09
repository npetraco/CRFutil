library(CRFutil)
library(MASS)
library(gRim)

lizards.info <- generate_brunner_lizards(typ = 2)
X.ct         <- lizards.info$lizards.contingency.table

# Specified graphs:
# Fit with dmod, Hojgaard gRim. Uses loglin
llm1 <- dmod(~1:2 + 1:3, data=X.ct)
#llm1$fitinfo$param

# Convert to crf object in standard parameterization
llm1.crf <- loglin2crf(llm1, standard.potentialsQ = T, plotQ = T)
dump.crf(llm1.crf)

# Original loglin() parameterization
llm1a.crf <- loglin2crf(llm1, standard.potentialsQ = F, plotQ = T)
dump.crf(llm1a.crf)

# Make sure original potentials from loglin give the same joint distribution
compute.full.distribution(llm1.crf)$joint.distribution
compute.full.distribution(llm1a.crf)$joint.distribution

