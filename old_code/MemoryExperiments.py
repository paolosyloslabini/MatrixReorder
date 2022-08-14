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

  
for m in [10,20,50,100]:
    for density in [0.01,0.05,0.1,0.2,0.5]:
        graph = csr.make_random_CSR(n,m,density)
        
        for sim_name in csr.similarities: 
            memories = []
            tau_mult = 1
            taus = np.linspace(0,1,20)
            for tau in taus:
                if  "hamming" in sim_name:
                        tau_mult = graph.M
                if "special" in sim_name:
                    tau_mult= 3
                
                tau *= tau_mult
                do_merge = True
                sim_measure = lambda x,y,z : csr.similarities[sim_name](x,y,z, tau);
                grouping = graph.group_by_sim(sim_measure,do_merge)
                in_density, block_rows, nz_blocks = graph.blocking(grouping, sim_measure)
                memory = block_rows + nz_blocks + graph.tot_nz()/in_density
                memories.append(memory)
            plt.plot(taus,memories, label = sim_name, marker = "x", markersize = 5)
        
        plt.legend()
        name = f"memory_comparison/memory_comparison_n{n}_m{m}_d_{density}"
        plt.yscale("log")
        plt.xlim(0,1)
        #plt.ylim()
        plt.xlabel("tau")
        plt.ylabel("memory")
        plt.title(name)
        plt.savefig(name + ".jpg", dpi = 1000)
        plt.show()