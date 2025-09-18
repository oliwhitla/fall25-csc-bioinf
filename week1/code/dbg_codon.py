import copy
from typing import Dict, Set, Optional, List

def reverse_complement(key: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complement[key[i]]
    return ''.join(key)


class Node:
    def __init__(self, kmer: str) -> None:
        
        # seed set with int type 
        self._children: set[int] = {-1}       # set[int]
        self._count: int = 0
        self.kmer: str = kmer
        self.visited: bool = False
        self.depth: int = 0
        self.max_depth_child: int = -1 # sentiel 

    def add_child(self, child_indx: int) -> None:
        self._children.add(child_indx)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = -1

    def get_count(self):
        return self._count

    def get_children(self):
        # only return if not sentiel 
        return [c for c in self._children if c != -1]

    def remove_children(self, target):
        # keeps children type consistent by filtering out the sentiel
        self._children = {c for c in self._children if c not in target and c != -1}


class DBG:
    def __init__(self, k: int, data_list):
        self.k = k
        # Nodes dict seeded with permanent sentinel -1
        self.nodes = {-1: Node("SENTINEL")}
        self.kmer2idx = {"SENTINEL": -1}
        self.kmer_count = 0

        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        # check data list
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                    self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def show_count_distribution(self):
        count = [0] * 30
        for idx, node in self.nodes.items():
            if idx == -1:    # skip sentiel
                continue
            c = node.get_count()
            # guard for out of bounds
            if 0 <= c < len(count):
                count[c] += 1
        print(count[0:10])



    def _add_node(self, kmer: str) -> int:
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str) -> None:
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child: int) -> int:
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx: int) -> list[int]:
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx):
        n = self.nodes[idx]
        if not n.visited:
            n.visited = True
            max_depth, max_child = 0, -1   # max depth says as int
            for child in self._get_sorted_children(idx):
                d = self._get_depth(child)
                if d > max_depth:
                    max_depth, max_child = d, child
            n.depth = max_depth + 1
            n.max_depth_child = max_child
        return n.depth

    def _reset(self) -> None:
        for idx, node in self.nodes.items():
            if idx == -1:   # skip sentinel
                continue
            node.reset()

    def _get_longest_path(self) -> list[int]:
        max_depth, max_idx = 0, -1  # sentinel instead of None
        for idx in self.nodes:
            if idx == -1:   # ✅ skip sentinel
                continue
            d = self._get_depth(idx)
            if d > max_depth:
                max_depth, max_idx = d, idx

        path = []
        while max_idx != -1:   # consistent int check
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: list[int]) -> None:
        path_set = set(path)
        for idx in path:
            if idx in self.nodes and idx != -1:
                del self.nodes[idx]
        for idx in list(self.nodes.keys()):
            if idx == -1:   # skip sentinel
                continue
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path: list[int]):
        if not path:
            return ""       # empty string instead of None
        concat = self.nodes[path[0]].kmer
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self):
        # reset params in nodes for getting longest path
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
