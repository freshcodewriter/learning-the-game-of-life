import os
from ete3 import Tree
from multiprocessing import Pool
from multiprocessing import cpu_count

DIR = "/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/data/perl-test/"
INPUT_DIR = DIR + "cell-names/"
OUT_DIR = "/home/derek/studyAboard/2019fall/ESE546/project/learning-the-game-of-life/data/perl-test/nw-tree/"
t = Tree('16div_200_targets_30_states_REF.txt', format=9)
file_list = []


def prune_tree(filename):
    print(filename)
    with open(filename) as f:
        node_list = [line.rstrip('\n') for line in f]
    file_base_name = os.path.basename(filename)
    t_prune = t.copy()
    t_prune.prune(node_list[0:3])
    t_prune.write(format=9, outfile=OUT_DIR + file_base_name + ".nw")


if __name__ == '__main__':
    for file in os.listdir(INPUT_DIR):
        filename = INPUT_DIR + file
        print(filename)
        file_list.append(filename)

    print(file_list)
    print('run in parallel')
    pool = Pool(processes=cpu_count()-1)
    results = pool.map(prune_tree, file_list)
    # pool.start()
    # pool.join()