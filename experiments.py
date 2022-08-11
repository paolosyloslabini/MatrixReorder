# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import csr
import similarities as sims
import cmat_reorderings as blk
import numpy as np


def evaluate_blocking(cmat, grouping, block_size):
    groups = np.unique(grouping)
    block_count = 0;
    total_block_height = 0;
    total_fill_in = 0;
    for group in groups:
        pattern = []
        group_size = 0

        #build the pattern;
        for row, g in enumerate(grouping):
            if g == group:
                pattern = sims.merge_patterns(pattern, cmat.pos[row])
                group_size += 1
        
                
        #count nz_blocks; add their size
        block_pattern = sims.get_block_pattern(pattern, block_size)
        block_count += len(block_pattern);
        total_block_height += len(block_pattern)*group_size
        
        #count nonzeros
        for row, g in enumerate(grouping):
            if g == group:
                #now count nonzeros
                blocked_row = sims.get_block_pattern(cmat.pos[row],block_size)
                total_fill_in += block_size*len(np.setdiff1d(block_pattern, blocked_row, True));
                total_fill_in += len(blocked_row)*block_size - len(cmat.pos[row])
        
    return block_count, total_block_height, total_fill_in




def run_experiments(cmat, blocking, dist_func, taus, block_size):
    densities = []
    avg_sizes = []
    
    for tau in taus:
        grouping = blocking(cmat, tau, dist_func)
        block_count, total_block_height, total_fill_in = evaluate_blocking(cmat, grouping, block_size)
        
        density = cmat.tot_nz()/(cmat.tot_nz() + total_fill_in)
        densities.append(density)
        avg_sizes.append(total_block_height/block_count)
    return densities,avg_sizes


cmat = csr.make_random_CSR(100,100,0.1);
taus = np.linspace(0, 1, 50);

block_size = 1;
dist_func = lambda p1,p2, size1: sims.JaccardGeneral(p1, p2, block_size, 1, 1)
results = run_experiments(cmat, blk.IterativeClusteringPattern,dist_func, taus, block_size)


