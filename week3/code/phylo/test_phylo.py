# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further

import time
import os
import numpy as np
# figure out how this actually works 
from tree import Tree, TreeNode
from upgma import upgma
from nj import neighbor_joining



def data_dir(subdir):
    """
    Return the absolute path to the 'data/<subdir>' folder
    located one level above 'code/'.
    """
    # Start from the folder this test file lives in
    base_dir = os.path.dirname(os.path.realpath(__file__))

    # Go up two levels: from code/phylo/ → project_root/
    project_root = os.path.abspath(os.path.join(base_dir, "..", ".."))

    # Join project root with data/<subdir>
    return os.path.join(project_root, "data", subdir)



def distances():
    path = os.path.join(data_dir("sequence"), "distances.txt")
    return np.loadtxt(path, dtype=int)



def upgma_newick():
    path = os.path.join(data_dir("sequence"), "newick_upgma.txt")
    with open(path, "r") as file:
        newick = file.read().strip()
    return newick



def tree(distances):
    return upgma(distances)




# # what is DendroUpGMA? 
def test_upgma(tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):

            tree_dist = tree.get_distance(i, j)
            ref_dist = ref_tree.get_distance(i, j)
            

            # manual approx with tolerance 1e-3
            if abs(tree_dist - ref_dist) >= 1e-3:
                raise AssertionError(
                    f"Distance mismatch at ({i}, {j}): {tree_dist} vs {ref_dist}"
                )
            
            tree_topo= tree.get_distance(i, j, topological=True)
            ref_topo= ref_tree.get_distance(i, j, topological=True)
            
            if tree_topo != ref_topo:
                raise AssertionError(
                    f"Topological distance mismatch at ({i}, {j}): {tree_topo} vs {ref_topo}"
                )




# # So the test checks whether your algorithmically built tree (neighbor_joining(dist)) 
# # matches the manually built one using your TreeNode structure.

def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = neighbor_joining(dist)

    assert test_tree == ref_tree






def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10


def main():
    distances_var = distances()
    newick = upgma_newick()
    tree_var = tree(distances_var)
    
    start_time = time.time()

    test_distances(tree_var)
    print("Distance Tree Passed\n")
    test_neighbor_joining()
    print("Neighbor Joining Passed\n")
    test_upgma(tree_var, newick)
    print("Upgma Passed")
    

    end_time = time.time()
    elapsed_ms = (end_time - start_time) * 1000
    
    print(f"python: {elapsed_ms:.0f}ms")



if __name__ == "__main__":
    main()




