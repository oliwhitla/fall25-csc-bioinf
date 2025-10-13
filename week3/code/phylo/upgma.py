# upgma.py
import numpy as np
from tree import Tree, TreeNode

def upgma(distances):
    """
    Perform hierarchical clustering using UPGMA 
    (Unweighted Pair Group Method with Arithmetic Mean)
    """

    distances = np.asarray(distances, dtype=float)

    # validation checks
    if distances.shape[0] != distances.shape[1]:
        raise ValueError("Distance matrix must be square")
    if not np.allclose(distances, distances.T):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances < 0).any():
        raise ValueError("Distances must be non-negative")

    n = distances.shape[0]
    nodes = [TreeNode(index=i) for i in range(n)]
    is_clustered = np.zeros(n, dtype=bool)
    cluster_size = np.ones(n, dtype=int)
    node_heights = np.zeros(n, dtype=float)

    D = distances.copy().astype(float)

    while True:
        # find minimum distance among active clusters
        dist_min = float("inf")
        i_min = j_min = -1

        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                if D[i, j] < dist_min:
                    dist_min = D[i, j]
                    i_min, j_min = i, j

        # Stop when no more pairs to merge
        if i_min == -1 or j_min == -1:
            break

        # merge clusters i_min and j_min
        height = dist_min / 2.0
        nodes[i_min] = TreeNode(
            [nodes[i_min], nodes[j_min]],
            [height - node_heights[i_min], height - node_heights[j_min]]
        )
        node_heights[i_min] = height

        # mark j_min as inactive
        nodes[j_min] = None
        is_clustered[j_min] = True

        # update distances using weighted mean
        total_size = cluster_size[i_min] + cluster_size[j_min]
        for k in range(n):
            if not is_clustered[k] and k != i_min:
                mean = (
                    D[i_min, k] * cluster_size[i_min] +
                    D[j_min, k] * cluster_size[j_min]
                ) / total_size
                D[i_min, k] = D[k, i_min] = mean

        cluster_size[i_min] = total_size

    # return the last non-None node as root
    for root in reversed(nodes):
        if root is not None:
            return Tree(root)

    raise ValueError("Failed to build UPGMA tree: input may be empty")
