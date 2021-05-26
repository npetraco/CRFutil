library(CRFutil)
library(gRbase)
library(Rgraphviz)
library(MASS)
library(gRim)

grf.eq <- ~A + B + A:B
AB <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)
AB$edge.pot[[1]][1,1] <- AB$edge.pot[[1]][2,2] <- 2 # ******* theta true 11,22 is theta-times more likely than 12 or 21

num.samps <- 500
samps     <- sample.exact(AB, num.samps)
mrf.sample.plot(samps)
#samps
head(samps)

X  <- samps                              # Raw Observed states
Xc <- xtabs(~., data=data.frame(X))                  # Contingency table
Xf <- as.data.frame(xtabs(~., data=data.frame(X)))   # Freq table of observed states
X
Xc
Xf
megi <- marginal.edge.loglin(samps,conf.level = 0.95)


llm1 <- loglm( ~ X1 + X2 + X1:X2, data = Xc)
llm1$param
megi$glm.theta.raw


# theta (energies) According to glm p-values of coefs:
megi$glm.theta.est

grf.eq <- ~A+B+A:B
AB <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)
AB$node.pot <- megi$glm.poi.rescaled.node.pot
AB$edge.pot[[1]] <- megi$glm.poi.rescaled.edge.pot
AB$node.pot #
AB$edge.pot

junk <- marginal.edge.bels(AB, node.names = c("A","B"), state.names = c("0","1"))
as.numeric(junk$`Bel(A,B)`)
tm <- cbind(c(1.1,2.1),c(1.2,2.2))
tm
as.numeric(tm)
tv <- c("1","2")
tv2 <- c("X1","X2")
paste0(tv2[1],"=",tv[1],",",tv2[2],"=",tv[1])
