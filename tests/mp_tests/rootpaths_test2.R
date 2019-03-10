library(CRFutil)

# Model:
grphf <- ~1:3 + 2:3 + 3:4 + 3:5 + 5:6 + 6:7

# Find all paths from root to leaves:
gplst <- get.root.paths(grphf,root.node = 6)
paths.to.serial.scheds
