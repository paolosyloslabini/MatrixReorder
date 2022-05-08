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

try:
    os.mkdir("similarity_curves")
except:
    0    
try:
    os.mkdir("block_structures")
except:
    0

class CSR: 
    def __init__(self):
        self.clean()
        
    def clean(self):
        self.N = 0
        self.M = 0
        self.nzcount = []
        self.pos = []
        
    def __str__(self):
        print("N:", self.N, "M:", self.M)
        print("nzcount: ", self.nzcount)
        print("pos: ", self.pos)
        return " "
        
    
    def tot_nz(self):
        tot = 0
        for count in self.nzcount:
            tot += count
        return tot
    
    def fill_from_edgelist(self,edgelist_file, delimiter = " "):
        #must be a sorted edgelist (gaps allowed)
        self.clean()
        with open(edgelist_file) as f:
            for line in f:
                linesplit = line.split(delimiter)
                inc = int(linesplit[0]);
                out = int(linesplit[1]);
                
                while inc > self.N - 1:
                    #create new node;
                    self.add_node()
                    
                self.add_edge(row = -1, col = out)
                
        return self.N;
    
    def fill_from_array(self,array):
        self.M = len(array[0])
        for row_idx, row in enumerate(array):
            self.add_node()
            for col_idx, elem in enumerate(row):
                if elem != 0.:
                    self.nzcount[row_idx] +=1;
                    self.pos[row_idx].append(col_idx)
                    
    
    def add_node(self):
        self.N += 1
        self.nzcount.append(0)
        self.pos.append([])
    
    def add_edge(self, row, col):
        self.nzcount[row] += 1;
        self.pos[row].append(col);
        
    def group_by_sim(self,sim_func, do_merge = True):
        
        group_name = -1;
        group_array = np.ones(self.N)*(-1)
        for row_idx in range(self.N):
            if group_array[row_idx] == -1:
                group_name += 1;
                group_size = 1
                group_array[row_idx] = group_name
                pattern = self.pos[row_idx]
                for other_row_idx in range(row_idx+1,self.N):
                    #print("comparing", row_idx, other_row_idx)
                    if group_array[other_row_idx] == -1:
                        other_pattern =  self.pos[other_row_idx]
                        merge = sim_func(group_size, pattern, other_pattern)
                        if merge:
                            group_size += 1;
                            group_array[other_row_idx] = group_name;
                            if do_merge:
                                pattern = merge_patterns(pattern, other_pattern)

        return group_array
    
    def blocking(self,grouping, sim, verbose = True):
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        
        current_group = 0;
        current_pattern = [];       
        
        new_idx = 0
        n_of_block_rows = 0
        nz_blocks = 0
        fake_nz = 0
        blocks_distribution = {}
        while new_idx < len(induced_row_order):
            old_idx = induced_row_order[new_idx]
            
            size = 1
            while new_idx < len(induced_row_order) and grouping[induced_row_order[new_idx]] == current_group:
                
                

                old_idx = induced_row_order[new_idx]
                row  = self.pos[old_idx]
                new_idx +=1
                size += 1
                
                old_zeros = len(np.setdiff1d(current_pattern, row, assume_unique = True))
                
                fake_nz += old_zeros + len(np.setdiff1d(row, current_pattern,assume_unique = True))
                                    
                current_pattern = np.union1d(current_pattern,row);
            
            if len(current_pattern) != 0:
                n_of_block_rows += 1
                nz_blocks += len(current_pattern)
                if size not in blocks_distribution:
                    blocks_distribution[size] = 1
                else:
                    blocks_distribution[size] += 1

            
            current_group = grouping[old_idx];
            current_pattern = [];
            
            
            density = self.tot_nz()/(self.tot_nz() + fake_nz)
            
            if (verbose):
                print("******************BLOCKING COMPLETED")
                print(f"BLOCK ROWS: {n_of_block_rows} of AVG. SIZE: {self.N/n_of_block_rows}")
                print(f"TRUE NONZEROS: {self.tot_nz()} FAKE NONZEROS : {fake_nz}, with AVG. in-block DENSITY: {density}")
                print("PRINTING BLOCK DISTRIBUTION: size -- blocks with that size")
            
                for num in sorted(blocks_distribution):
                    print(f"{num} --> {blocks_distribution[num]}")
            
        return density, n_of_block_rows, nz_blocks

    
    
    def blocking_show(self,grouping, sim, filename = "blocking_structure_example.txt"):
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        
        current_group = 0;
        current_pattern = [];
           
        
        with open(filename, "w") as outfile:
            new_idx = 0
            n_of_block_rows = 0
            nz_blocks = 0
            fake_nz = 0
            blocks_distribution = {}
            while new_idx < len(induced_row_order):
                old_idx = induced_row_order[new_idx]
                
                size = 1
                while new_idx < len(induced_row_order) and grouping[induced_row_order[new_idx]] == current_group:
                    
                    

                    old_idx = induced_row_order[new_idx]
                    row  = self.pos[old_idx]
                    new_idx +=1
                    size += 1
                    
                    old_zeros = len(np.setdiff1d(current_pattern, row, assume_unique = True))
                    
                    fake_nz += old_zeros + len(np.setdiff1d(row, current_pattern,assume_unique = True))
                    
                    printrow = [0 for _ in range(old_zeros)] + [1 for _ in range(len(row))]
                    
                    print(printrow, file = outfile)
        
                    current_pattern = np.union1d(current_pattern,row);
                
                if len(current_pattern) != 0:
                    n_of_block_rows += 1
                    nz_blocks += len(current_pattern)
                    if size not in blocks_distribution:
                        blocks_distribution[size] = 1
                    else:
                        blocks_distribution[size] += 1

                
                outfile.write("\n")
                current_group = grouping[old_idx];
                current_pattern = [];
            
            
            density = self.tot_nz()/(self.tot_nz() + fake_nz)
            print("******************BLOCKING COMPLETED", file = outfile)
            print(f"BLOCK ROWS: {n_of_block_rows} of AVG. SIZE: {self.N/n_of_block_rows}", file = outfile)
            print(f"TRUE NONZEROS: {self.tot_nz()} FAKE NONZEROS : {fake_nz}, with AVG. in-block DENSITY: {density}", file = outfile)
            print("PRINTING BLOCK DISTRIBUTION: size -- blocks with that size", file = outfile)
            
            for num in sorted(blocks_distribution):
                print(f"{num} --> {blocks_distribution[num]}", file = outfile)
            
            
        return density, n_of_block_rows, nz_blocks
    
    def print_blocking(self, grouping):
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        current_group = 0;
        
        for new_idx,old_idx in enumerate(induced_row_order):
            if grouping[old_idx] != current_group:
                print("\n")
                current_group = grouping[old_idx];
            print(self.pos[old_idx])
            
            
#****************************************************************

def print_pattern(pattern,m):
     #tofix
    to_print = ""
    last = 0
    for elem in pattern:
        for i in range(elem - last):
            to_print += "0 "
        last = elem
        to_print += "1 "
    for elem in range(m - pattern[-1]):
        to_print += "0 "

    print(to_print)
def merge_patterns(p1,p2):
    #can be made faster since they are sorted
    return np.union1d(p1, p2)


def weighted_sim(size1,p1,p2, tau = 0.7, use_size = False, relative_val = True, cosine = False):
    if (len(p1) == 0) and (len(p2) == 0):
        return True

    if (len(p1) == 0) or (len(p2) == 0):
        return False
    i = 0
    j = 0 
    unsim = 0;
    tot = 0;
    while i < len(p1) and j < len(p2):
        if p1[i] < p2[j]:
            i += 1
            tot += 1
            unsim += 1
        elif p1[i] > p2[j]:
            j += 1
            tot += 1
            if (use_size):
                unsim += size1
            else:
                unsim += 1
        else:
            tot += 1
            i += 1
            j += 1
        
    while i < len(p1):
        i += 1
        tot += 1
        unsim += 1
    while j < len(p2):
        j += 1
        tot += 1
        if use_size:
            unsim += size1
        else:
            unsim += 1
    if relative_val:   
        unsim = unsim*1./tot
    if cosine:
        unsim = unsim*1./(len(p1)*len(p2))
    return (unsim < tau)

    

def make_random_CSR(n,m, density):
    graph = CSR();
    k = int(density*n*m)
    
    nz_pos = np.random.choice(n*m, k)
    mat = np.zeros(n*m);
    mat[nz_pos] = 1
    mat = mat.reshape((n,m))
    graph.fill_from_array(mat)

    return graph

def run_experiments(graph, sim, taus, do_merge = True):
    fill_ins = []
    avg_sizes = []
    
    for tau in taus:
        sim_measure = lambda x,y,z : sim(x,y,z, tau);
        grouping = graph.group_by_sim(sim_measure, do_merge = do_merge)
        fake_nz, number_of_groups = graph.blocking(grouping)
        fill_ins.append(fake_nz)
        avg_sizes.append(graph.N/number_of_groups)
    return fill_ins,avg_sizes

def plot_comparison(n,m,density, folder = "similarity_curves"):
    plt.legend()
    name = f"{folder}/similarity_comparison_n{n}_m{m}_d_{density}"
    plt.yscale("log")
    plt.xlim(0.1,1)
    plt.ylim(1,n)
    plt.xlabel("density")
    plt.ylabel("avg_group_size")
    plt.title(name)
    plt.savefig(name + ".jpg", dpi = 1000)
    plt.show()
    
    
jaccard = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = True)
jaccard_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = True)

hamming = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = False)
hamming_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = False)

cosine = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = False, cosine = True)
cosine_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = False, cosine = True)

similarities = {
    "jaccard" : jaccard,
    "jaccard_special" : jaccard_special,
    "hamming": hamming,
    "hamming_special": hamming_special,
    "cosine" : cosine,
    "cosine_special": cosine_special
    }