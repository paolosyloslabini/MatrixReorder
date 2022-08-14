# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:14:55 2022

@author: Paolo
"""


import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import csr
from scipy.stats import binom


def prob_of_better(n, row1, row2, density):
    diff = len(np.setdiff1d(row2, row1))
    return binom.cdf(diff - 1,n - len(row1), density)

def merge_prob(graph, candidates, pattern, target_size, row):
    if len(pattern) == 0 and len(row) == 0:
        return True
    if len(row) == 0 or len(pattern) == 0:
        return False
    prob = prob_of_better(graph.M, pattern, row , graph.density())
    merge =  int(prob*candidates) < target_size*2
    #print(pattern, row, candidates, prob, merge)
    return merge    
    
def fixed_size_blocking(graph, block_size):
    group_name = -1;
    group_array = np.ones(graph.N)*(-1)
    
    for row_idx in range(graph.N):
        if group_array[row_idx] == -1:
            group_name += 1;
            group_size = 1
            group_array[row_idx] = group_name
            pattern = graph.pos[row_idx]
            size_reached = False
            for other_row_idx in range(row_idx+1,graph.N):
                if size_reached:
                    break
                if group_array[other_row_idx] == -1:
                    other_pattern =  graph.pos[other_row_idx]
                    candidates = len([x for x in group_array[other_row_idx + 1:] if x == -1])
                    merge = merge_prob(graph, candidates, pattern, block_size - group_size, other_pattern)
                    if merge:
                        group_size += 1;
                        candidates = -1;
                        group_array[other_row_idx] = group_name;
                        pattern = csr.merge_patterns(pattern, other_pattern)
                        if group_size == block_size:
                            size_reached = True
    return group_array
                        




n = 1024
m = 100
density = 0.1

graph = csr.make_random_CSR(n, m, density)
grouping = fixed_size_blocking(graph, 10)
#graph.print_blocking(grouping)
graph.blocking(grouping, verbose = True)