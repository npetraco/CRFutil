library(CRFutil)

# Model:
#grf <- ~1:2# Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX

# Convert to pw factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)

# Schedules:
schs <- get.root.paths(pwfg, root.node = 4, serial.schedsQ = T)
schs$forward
schs$backward

# Initialize a storage list to hold messages:
msg.cont.info <- init.message.storage(schs)
msg.bxs       <- msg.cont.info$message.container
msg.cont.info$message.names.mat
msg.cont.info$message.names.mat[,1]
msg.cont.info$message.names.mat[,2]
msg.cont.info$message.names.mat[,3]


msg.nms       <- names(msg.bxs)
msg.nms[1]
msg.nms[2]

# message type test:
msg.nms[2]
message.type(msg.nms[2])
message.type("1")
message.type("f1")

# Neighbors except for, nex, test
nex(pwfg, "f34", "4") # find all neighbors of f34 except 4
nex(pwfg, "4", "f34") # find all neighbors of 4 except f34
nex(pwfg, "f45", "5")

# get messages test
msg.nms[2]
strsplit(msg.nms[2], split = ".", fixed = T)[[1]]
#

paste0(nex(pwfg, "4", "f34"), ".", "4")

get.incoming.messages(
  start.node="4",
  neighborx.nodes=nex(pwfg, "4", "f34"), ".", "4"),

  message.list=msg.bxs)

s# f34 -> 4
adj(pwfg, "f34")
# 4 -> f24
