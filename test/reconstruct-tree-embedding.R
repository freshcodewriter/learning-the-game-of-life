library(phangorn)
library(castor)
DATA <- read.csv(file="/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/data/dataset-3/simulations/07div_200_targets_30_states_0378rep_FAST_0077seed.vec", header=FALSE, sep=",")

rownames(DATA) = 1:128
colnames(DATA) = 1:200

dm_2 = matrix(0, 128, 128)
for (i in c(1: 128))
{
  row_matrix = matrix(rep(DATA[i,], nrow=128, byrow = TRUE))
  bool_matrix = DATA - row_matrix != 0
  dist_col = rowSums(bool_matrix)
  dm_2[, i] = dist_col
}
my_dim <- as.dist(dm_2)

tree <- upgma(my_dim)
tree$edge.length <- rep(1, 254)
plot(tree, font = 1, cex = 0.5)