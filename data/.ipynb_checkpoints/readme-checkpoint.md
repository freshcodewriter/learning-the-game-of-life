
# Learning The Game of Life Dataset

## Directory Structure

```bash
.
├── dataset-1
│  ├── ref-trees
│  │  └── 16div_200_targets_30_states_REF.nw
│  └── simulations
│     ├── 16div_200_targets_30_states_0001rep_FAST.fas
│     ├── ...
│     └── 16div_200_targets_30_states_0010rep_FAST_freqs.txt
├── dataset-2
│  ├── ref-trees
│  │  └── 10div_200_targets_30_states_REF.nw
│  └── simulations
│     ├── 10div_200_targets_30_states_0001rep_FAST.fas
│     ├── ...
│     └── 10div_200_targets_30_states_0010rep_FAST_freqs.txt
└── scripts
   └── treeGeneration.m
```

## Dataset desciption

### dataset-1

Dataset-1 contains 1 linage tree with 16 cell divisions, i.e. there are 65536 cells. Based on this linage tree, there are 10 difference possible cell sequencing data simulated. Generate one sequencing data cost < 3 mins so far.

## TODO

- Subset the linage tree to 1000 cells. Based on the paper, ideally 1000 cell should be able to reconstruct linage tree.
- Improve the matlab script for parallel simulating.
