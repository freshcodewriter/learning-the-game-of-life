DATA <- read.csv(file="/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/calculated_pairwise_matrices3/pairwise33.csv", header=TRUE, sep=",")
DATA <- DATA[,-1]
dim(DATA)
rownames(DATA) = 1:128
colnames(DATA) = 1:128
my_dim <- as.dist(DATA)
tree <- upgma(my_dim)
tree$edge.length <- rep(1, 254)
plot(tree, font = 1, cex = 0.5)

# ref_tree = read_tree(file='/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/data/dataset-3/ref-trees/07div_200_targets_30_states_REF_trim_n.nw',
#                      include_edge_lengths=FALSE,
#                      check_label_uniqueness=TRUE)
# ref_tree$tip.label <- c(128:1)
# ref_tree$edge.length <- rep(1, 254)
# 
# 
# res <- comparePhylo(tree, ref_tree)
# x <- strsplit(res$messages, '\t')[[8]]
# x <- as.character(x)
# x <- strsplit(x, ' ')[[1]][1]
# x <- as.numeric(x)
# y <- 127 - x
# print(y)