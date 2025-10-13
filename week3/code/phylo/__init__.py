"""
Python port of Biotites `biotite.sequence.phylo` module.

This simplified version reimplements selected functionality from the
original Cython sources (`tree.pyx`, `upgma.pyx`, `nj.pyx`) in pure
Python for compatibility with Codon.

Only the classes and functions required for the Week 3 deliverable
are included:
    - Tree
    - TreeNode
    - upgma()
    - neighbor_joining()

Original Biotite source: https://github.com/biotite-dev/biotite
License: 3-Clause BSD
"""


from .tree import Tree, TreeNode
from .upgma import upgma
from .nj import neighbor_joining

__all__ = ["Tree", "TreeNode", "upgma", "neighbor_joining"]
