#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly=TRUE)

print(args)

tree = ape::read.tree(args[1])

tips = read.csv(args[2])$sequence_name

#tree$tip.label = sapply(strsplit(tree$tip.label,"/"), `[`, 2)

print("tips not in:")
print(tips[!(tips %in% tree$tip.label)])
print(sum(!(tips %in% tree$tip.label)))

tips = tips[tips %in% tree$tip.label]

tips_not_in = tree$tip.label[!(tree$tip.label %in% tips)]

subtree = drop.tip(tree,tips_not_in)

print("got subtree")

subtree = multi2di(subtree)

write.tree(subtree, file=args[3])
