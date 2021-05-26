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

samp.info  <- sim.field.random(adj, n.states, 100)
true.model <- samp.info$model
true.model$node.pot
true.model$edge.pot


ft.ifo <- mrf.standard.fit(samp.info$samples, grphf, n.states, mrf.exact.nll, infer.exact)
ft.ifo <- mrf.standard.fit(samp.info$samples, grphf, n.states, mrf.junction.nll, infer.junction)

pot.info <- make.gRbase.potentials(ft.ifo$fit.model, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

bel.info <- make.gRbase.beliefs(ft.ifo$inference.info, node.names = gp@nodes, edge.mat = ft.ifo$fit.model$edges, state.nmes=NULL)
bel.info$node.beliefs[[1]]
ft.ifo$inference.info$node.bel[1,]
bel.info$node.beliefs[[2]]
ft.ifo$inference.info$node.bel[2,]
bel.info$node.beliefs[[3]]
ft.ifo$inference.info$node.bel[3,]
bel.info$node.beliefs[[4]]
ft.ifo$inference.info$node.bel[4,]

bel.info$edge.beliefs[[1]]
ft.ifo$inference.info$edge.bel[[1]]
bel.info$edge.beliefs[[2]]
ft.ifo$inference.info$edge.bel[[2]]
bel.info$edge.beliefs[[3]]
ft.ifo$inference.info$edge.bel[[3]]
