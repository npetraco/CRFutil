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

# Grab messages of requisite neighboring nodes needed to construct the outgoing message
# This task needs msg.cont.info and the factor graph
msg.num      <- 16
msg.sch[msg.num,]
msg.sp       <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)

# st.nex.end   <- nex(pwfg, msg.sp[1], msg.sp[2])       # neighbors of start node except end node. Needs factor graph
# msg.incoming <- paste0(nex(pwfg, msg.sp[1], msg.sp[2]) , ".", msg.sp[1]) # incoming messages to start node
# msg.incoming
#
# msg.incoming.idxs <- sapply(1:length(msg.incoming), function(xx){which(msg.sch[,4] == msg.incoming[xx])})
# msg.incoming.list <- msg.bxs[msg.incoming.idxs]
# msg.incoming.list  # For test function, these should be assigned to "Done", so check for that. Can we use the "Done" sequence to further optimize the schedule??
# names(msg.incoming.list) # requisite incoming messages from neighboring nodes needed to construct the outgoing message
# msg.incoming

