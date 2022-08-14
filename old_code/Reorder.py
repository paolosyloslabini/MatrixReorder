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
    

n = 1000;
m = 20;
density = 0.2;
marker = itertools.cycle(('*',',', '+', '.', 'o', '*')) 

  
for m in [10,20,50,100]:
    for density in [0.01,0.05,0.1,0.2,0.5]:
        graph = csr.make_random_CSR(n,m,density)
        
        
     
        for sim_name in csr.similarities: 
            densities = []
            avg_sizes = []   
            for tau in np.linspace(0,1,20):
                if  "hamming" in sim_name:
                        tau *= graph.M
                if "special" in sim_name:
                    tau *= 3
                do_merge = True
                sim_measure = lambda x,y,z : csr.similarities[sim_name](x,y,z, tau);
                grouping = graph.group_by_sim(sim_measure,do_merge)
                in_density, number_of_blocks = graph.blocking_show(grouping, sim_measure, filename = f"block_structures/block_structure_n{n}_m{m}_d{density}_{sim_name}_t{tau}.txt")
                densities.append(in_density)
                avg_sizes.append(graph.N/number_of_blocks)
            plt.plot(densities,avg_sizes, label = sim_name, marker = next(marker))
        csr.plot_comparison(n, m, density)