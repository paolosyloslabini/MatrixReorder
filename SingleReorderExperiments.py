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
from matplotlib.lines import Line2D


n = 1000;
m = 20;
density = 0.2;

taus = np.linspace(0,1,20)
repetitions = 10
for m in [10,20,50,100]:
    for density in [0.01,0.05,0.1,0.2,0.5]:
        
        for sim_name in csr.similarities:
            
            tau_mult = 1
            if  "hamming" in sim_name:
                    tau_mult = m
            if "special" in sim_name:
                tau_mult= 3
            
            marker = itertools.cycle(Line2D.markers)
            for tau in taus:
                tau *= tau_mult
                
                densities = []
                sizes = []
                
                for rep in range(repetitions):
                    graph = csr.make_random_CSR(n,m,density)
                    do_merge = True
                    sim_measure = lambda x,y,z : csr.similarities[sim_name](x,y,z, tau);
                    grouping = graph.group_by_sim(sim_measure,do_merge)
                    in_density, number_of_blocks = graph.blocking_show(grouping, sim_measure, filename = f"block_structures/block_structure_n{n}_m{m}_d{density}_{sim_name}_t{tau}.txt")
                    densities.append(in_density)
                    sizes.append(graph.N/number_of_blocks)
                
                
                plt.scatter(densities,sizes, label = tau, marker = next(marker))

            plt.legend()
            name = f"multiple_mat_exp/repeated_{sim_name}_n{n}_m{m}_d_{density}"
            plt.yscale("log")
            plt.xlim(0,1)
            #plt.ylim()
            plt.xlabel("tau")
            plt.ylabel("memory")
            plt.title(name)
            plt.savefig(name + ".jpg", dpi = 1000)
            plt.show()