jmat <- t(array(c("A","B","C","D","E","FF","G","H","II","J","K","L","M","N","O","P") ,c(4,4)))
jmat

geq <- make.lattice(num.rows = 4, num.cols = 4, cross.linksQ = T, plotQ = T)
gp <- ug(geq, result = "graph")
plot(gp)
dev.off()

geq <- make.lattice(num.rows = 4, num.cols = 4, cross.linksQ = F, plotQ = T)
gp <- ug(geq, result = "graph")
plot(gp)
dev.off()

geq <- make.lattice(num.rows = 4, num.cols = 4, cross.linksQ = T, plotQ=T, node.names.vec = c("A","B","C","D","E","FF","G","H","II","J","K","L","M","N","O","P"))
gp <- ug(geq, result = "graph")
plot(gp)
dev.off()

geq <- make.lattice(num.rows = 4, num.cols = 4, cross.linksQ = F, plotQ=T, node.names.vec = c("A","B","C","D","E","FF","G","H","II","J","K","L","M","N","O","P"))
gp <- ug(geq, result = "graph")
plot(gp)
dev.off()


geq <- make.lattice(num.rows = 32, num.cols = 32, cross.linksQ = T)
gp <- ug(geq, result = "graph")
dev.off()
plot(gp)
dev.off()


geq <- make.lattice(num.rows = 4, num.cols = 4, cross.linksQ = F)
adj        <- ug(geq, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)

fit <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
