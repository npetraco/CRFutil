library(CRFutil)

# Model:
#grf <- ~1:2# Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX
grf <- ~1:2 + 1:4 + 2:5 + 2:6 + 3:7 + 4:7 + 4:8 + 6:9 + 7:10 + 11:10 + 10:14 + 10:13 + 10:12 + 12:17 + 12:16 + 12:15

# Convert to pw factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = 10, serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch

# Grab messages of requisite neighboring nodes needed to construct the outgoing message
# This task needs msg.cont.info and the factor graph
msg.num      <- 52
msg.sch[msg.num,]
msg.sp       <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)

