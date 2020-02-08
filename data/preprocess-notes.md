# Preprocess Workflow

## treeGeneration.m 

generate tree in /Ref_Trees.
generate cell sequence in /simulations.
Simulated result could be used as dataset 1.

## Nexux_fromSims_AlleleCount.pl

Sample 1000 leaves from the output of matlab

## extract_cell-id.sh

extract the selected cell ids for prune reference tree.

## prune-test.py

generated sub-tree only containing the selected cells as its leaves.
