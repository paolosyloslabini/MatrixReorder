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

try:
    os.mkdir("similarity_curves")
except:
    0    
try:
    os.mkdir("block_structures")
except:
    0
    

n = 200;
m = 10;
  
for m in [20,50,100]:
    #for density in [0.05,0.1,0.2,0.5]:
    for density in [0.1,0.2]:
        graph = csr.make_random_CSR(n,m,density)
        

        for sim_name in ["special", "jaccard"]: 
            densities = []
            avg_sizes = []   
            for tau in np.linspace(0,1,20):
                do_merge = True
                sim_measure = lambda n1,v1,n2,v2: csr.similarities[sim_name](n1,v1,n2,v2) <= tau;
                grouping = graph.group_by_sim(sim_measure,do_merge)
                try:
                    os.mkdir(f"block_structures/{sim_name}")
                except:
                    None
                in_density, total_block_height, nz_blocks = graph.blocking_show(grouping, filename = f"block_structures/{sim_name}/block_structure_n{n}_m{m}_d{density}_{sim_name}_t{tau}.txt")
                densities.append(in_density)
                avg_sizes.append(total_block_height/nz_blocks)
            plt.plot(densities,avg_sizes, label = sim_name, marker = "x", markersize = 5)
        
        densities = []
        avg_sizes = []
           
        """
        for block_size in [2,4,8,16,32,64]:
            sim_measure = lambda y,z : csr.general_sim(1,y,z,use_size = False, relative_val = True)
            grouping = graph.hierchical_blocking_link(block_size, sim_measure)
            sim_name = "hierarchical"
            try:
                os.mkdir(f"block_structures/{sim_name}")
            except:
                None
            in_density, total_block_height, nz_blocks = graph.blocking_show(grouping, filename = f"block_structures/{sim_name}/block_structure_n{n}_m{m}_d{density}_{sim_name}_b{block_size}.txt")
            densities.append(in_density)
            avg_sizes.append(total_block_height/nz_blocks)
        plt.plot(densities,avg_sizes, label = sim_name, marker = "x", color = "r", markersize = 5)
        """

        """
        densities = []
        avg_sizes = []
        
        tot_zeros = graph.N*graph.M - graph.tot_nz()
        for htau in np.linspace(0,1,10):
            grouping = graph.hierchical_blocking_nofix(csr.gSpecial, htau)
            sim_name = "hierarchical_true"
            try:
                os.mkdir(f"block_structures/{sim_name}")
            except:
                None
            in_density, total_block_height, nz_blocks = graph.blocking_show(grouping, filename = f"block_structures/{sim_name}/block_structure_n{n}_m{m}_d{density}_{sim_name}_t{htau}.txt")
            densities.append(in_density)
            avg_sizes.append(total_block_height/nz_blocks)
        plt.plot(densities,avg_sizes, label = sim_name, marker = "x", color = "g", markersize = 5)
        """        
        
        """
        densities = []
        avg_sizes = []
        
        for merge_name, merge_crit in (("fixed_size",csr.merge_prob),):
            for tau in [2,4,8,16,32,64]:
                grouping = graph.fixed_size_blocking(tau, merge_crit)
            try:
                os.mkdir(f"block_structures/{sim_name}")
            except:
                None           
            in_density, total_block_height, nz_blocks = graph.blocking_show(grouping, filename = f"block_structures/{sim_name}/block_structure_n{n}_m{m}_d{density}_{sim_name}_b{block_size}.txt")
            densities.append(in_density)
            avg_sizes.append(total_block_height/nz_blocks)
        plt.plot(densities,avg_sizes, label = merge_name, marker = "x", markersize = 5)
        """
            
        csr.plot_comparison(n, m, density)