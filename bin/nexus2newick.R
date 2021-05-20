library(ape)

args = commandArgs(trailing=TRUE)

t <- read.nexus(args[1])

write.tree(t, args[2])
