ref_tree = read_tree(file='/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/data/dataset-3/ref-trees/07div_200_targets_30_states_REF_trim_n.nw',
                     include_edge_lengths=FALSE,
                     check_label_uniqueness=TRUE)
ref_tree$tip.label <- c(128:1)
ref_tree$edge.length <- rep(1, 254)
plot(ref_tree, font = 1, cex = 0.5)