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
density = 0.2;

taus = np.linspace(0,1,20)
repetitions = 20
for m in [10,20,50,100]:
    for density in [0.01,0.05,0.1,0.2,0.5]:
        
        for sim_name in csr.similarities:
            
            tau_mult = 1
            if  "hamming" in sim_name:
                    tau_mult = m
            if "special" in sim_name:
                tau_mult= 3
                
                
                
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
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
                    in_density, block_rows, nz_blocks = graph.blocking(grouping, sim_measure, verbose = False)
                    densities.append(in_density)
                    sizes.append(graph.N/block_rows)
                
                
                sc = plt.scatter(densities,sizes, c = [tau]*20, marker = next(marker), s = 2)
                plt.clim(0, tau_mult)


            name = f"repeated_exp/{sim_name}_n{n}_m{m}_d_{density}"
            plt.colorbar()
            plt.yscale("log")
            plt.xlim(0,1)
            plt.ylim(1,n)
            plt.xlabel("density")
            plt.ylabel("avg_size")
            plt.title(name)
            plt.savefig(name + ".jpg", dpi = 1000)
            plt.show()