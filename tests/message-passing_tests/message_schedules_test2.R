library(CRFutil)

# Model:
#grf <- ~1:2# Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX

# Convert to pw factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = 4, serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch


# message type test:
colnames(msg.sch)
msg.sch[,4]
message.type(msg.sch[1,4])
message.type(msg.sch[7,4])
message.type("1")
message.type("f1")

# Neighbors except for end node test:
adj(object = pwfg, index = "f34")
nex(pwfg, "f34", "4") # find all neighbors of f34 except 4

adj(object = pwfg, index = "4")
nex(pwfg, "4", "f34") # find all neighbors of 4 except f34

adj(object = pwfg, index = "f45")
nex(pwfg, "f45", "5") # find all neighbors of f34 except 4


# get messages test
stj <- "4"
spj <- "f34"
reqj <- paste0(nex(pwfg, stj, spj), ".", stj)

msg.sch[,4]
reqj[2]
which(msg.sch[,4] == reqj[2])

sapply(1:length(reqj), function(xx){which(msg.sch[,4] == reqj[xx])})

adj(object = pwfg, index = "f1")
nex(pwfg, "f1", "1") # find all neighbors of f1 except 1
length(nex(pwfg, "f1", "1")) # for leafs
