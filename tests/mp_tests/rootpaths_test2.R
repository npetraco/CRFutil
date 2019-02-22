library(CRFutil)

# Model:
grphf <- ~1:3 + 2:3 + 3:4 + 3:5 + 5:6 + 6:7

# Find all paths from root to leaves:
get.root.paths(grphf,root.node = 5)
