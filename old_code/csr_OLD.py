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
from scipy.stats import binom


def prob_of_better(n, row1, row2, density):
    diff = len(np.setdiff1d(row2, row1))
    return binom.cdf(diff - 1,n - len(row1), density)

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
    
    def density(self):
        return self.tot_nz()*1./(self.N*self.M)
    
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
                        merge = sim_func(pattern, group_size, other_pattern, 1)
                        if merge:
                            group_size += 1;
                            group_array[other_row_idx] = group_name;
                            if do_merge:
                                pattern = merge_patterns(pattern, other_pattern)

        return group_array
    
    
    
    def hierchical_blocking_link(self, target_block_size, similarity, current_block_size = 1, grouping = None):
        
        if current_block_size == 1:
            grouping = range(self.N)
        
        new_group = np.ones(self.N)*(-1)
        pointers = np.ones(self.N//current_block_size)*(-1)
        
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        
        #print(f"fixed size blocking, current size {current_block_size}, target {target_block_size}") 
        #print(f"CURRENT GROUPING: {grouping}")
    
        def is_assigned(group_idx):
            true_idx = group_idx*current_block_size
            old_idx = induced_row_order[true_idx]
            #print(f"checking {i}, old idx {old_idx}")
            return new_group[old_idx] != -1
        
        def assign_group(group_idx, group_name):
            true_idx = group_idx*current_block_size
            #print(f"assigning {i} to {group_name}")
            for idx in range(true_idx, true_idx + current_block_size):
                old_idx = induced_row_order[idx]
                new_group[old_idx] = group_name
                #print(f"***** {idx}, old_idx {old_idx}")
                
        def get_pattern(group_idx):
            
            true_idx = group_idx*current_block_size
            pattern = self.pos[induced_row_order[true_idx]]
            for idx in range(true_idx, true_idx+current_block_size):
                row = self.pos[induced_row_order[true_idx]]
                pattern = np.union1d(pattern,row);
            return pattern

        def find_unassigned():
            i = 0;
            while i < self.N//current_block_size and is_assigned(i):
                i += 1;
            if i == self.N//current_block_size: 
                return -1
            else:
                return i
            
        def find_min_idx(group_idx):
            min_idx = int(pointers[group_idx])
            if pointers[group_idx] == -1 or is_assigned(min_idx):
                min_sim = float('inf')
                min_idx = -1
                #print("recalculating min for {group_idx}")
                for j in range(self.N//current_block_size):
                    
                    if not is_assigned(j) and j != group_idx:
                        other_pattern = get_pattern(j)
                        tmp_sim = similarity(pattern, other_pattern)
                        if tmp_sim < min_sim:
                            min_sim = tmp_sim
                            min_idx = j
            return min_idx
            
            


        i = 0;
        last_visited = -1
        current_group_name = 0
        not_done = True
        
        while (not_done):
            
            pattern = get_pattern(i)
        
        
            #find nearest neighbour. recaculate if necessary
            min_idx = find_min_idx(i)
            pointers[i] = min_idx
            
            #if the last two are each other's NN, merge them
            if min_idx == last_visited:
                assign_group(i,current_group_name)
                assign_group(min_idx, current_group_name)
                current_group_name += 1
                last_visited = -1
                i = find_unassigned()
                if i == -1:
                    not_done = False
            else:
                last_visited = i
                i = min_idx
                
        
        current_block_size *= 2;
        if current_block_size >= target_block_size:
            return new_group
        else:
            return self.hierchical_blocking(target_block_size, similarity, current_block_size, new_group)

    def hierchical_blocking_nofix(self, similarity, tau = 0.8):
        
        class Group():
             def __init__(self, graph, row_idx):
                 self.rows = [row_idx,]
                 self.ID = row_idx
                 self.graph = graph
                 self.pattern = graph.pos[row_idx]
                 self.closest = -1
                 self.closest_val = float('inf')
             def add(self,group):

                 self.rows = group.rows + self.rows
                 self.pattern = np.union1d(self.pattern, group.pattern)


             def __eq__(self, other):
                 if isinstance(other, Group):
                     return self.ID == other.ID
                 return False
             def update_closest(self, sender):
                 new_val = self.distance(sender)
                 if new_val < self.closest_val:
                     self.closest = sender
                     self.closest_val = new_val
                     return True
                 return False
             
             def distance(self, other):

                 p1 = self.pattern
                 p2 = other.pattern
                 if len(p1) == 0 and len(p2) == 0:
                     return 0
                 s1 = len(self.rows)
                 s2 = len(other.rows)
#                 if s1 > 10 or s2 > 10:
#                     return 1;
                 d = similarity(p1,s1,p2,s2)
                 return d
                 
             def find_closest(self, group_list):
                 min_sim = float('inf')
                 best_group = None
                 for g in group_list:
                     if g.ID != self.ID:
                         tmp_sim = self.distance(g)
                         if tmp_sim < min_sim:
                             min_sim = tmp_sim
                             best_group = g
                 self.closest = best_group
                 self.closest_val = min_sim
                 if best_group == None:
                     print(f"DEBUG: grouplist size {len(group_list)}, group {self.ID} with pattern {self.pattern}, closest val {min_sim}")
                 
        groups = [Group(self, row) for row in range(self.N)]
        for g in groups:
            g.find_closest(groups)
        
        
        min_val = -1
        while(len(groups) > 2):
            #find closest elements
            min_val = float('inf')
            min_idx = -1
            for idx, g in enumerate(groups):
                if g.closest_val < min_val:
                    min_val = g.closest_val
                    min_idx = idx
                    print(tau,min_val)
            if min_val > tau:
                break;
            little = groups.pop(min_idx)
            big = little.closest
            big.add(little)
            for g in groups:
                if g.ID == big.ID:
                    g.find_closest(groups)
                    continue
                if g.closest.ID == big.ID or g.closest.ID == little.ID:
                    #if the new group is not the closest, run the closest search
                    if not g.update_closest(big):
                        g.find_closest(groups)
                else: 
                    g.update_closest(big)
                    
            
        grouping = np.ones(self.N)
        current = 0
        for g in groups:
            for row in g.rows:
                grouping[row] = current
            current += 1
        return grouping
        
    
    def hierchical_blocking(self, target_block_size, similarity, current_block_size = 1, grouping = None):
        
        if current_block_size == 1:
            grouping = range(self.N)
        
        new_group = np.ones(self.N)*(-1)
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        
        #print(f"fixed size blocking, current size {current_block_size}, target {target_block_size}") 
        #print(f"CURRENT GROUPING: {grouping}")
    
    
        def is_not_assigned(i):
            old_idx = induced_row_order[i]
            #print(f"checking {i}, old idx {old_idx}")
            return new_group[old_idx] == -1
        
        def assign_group(i, group_name):
            #print(f"assigning {i} to {group_name}")
            for idx in range(i, i + current_block_size):
                old_idx = induced_row_order[idx]
                new_group[old_idx] = group_name
                #print(f"***** {idx}, old_idx {old_idx}")
                
        def get_pattern(idx):
            #print(f"pattern:, {idx}, {len(induced_row_order)}")

            #print(f"getting pattern for row {idx}")
            pattern = self.pos[induced_row_order[idx]]
            for idx in range(idx, idx+current_block_size):
                row = self.pos[induced_row_order[idx]]
                pattern = np.union1d(pattern,row);
            return pattern

        
        current_group_name = 0
        for i in range(self.N//current_block_size):
            group_idx = i*current_block_size
            #print(f"EVALUATING row {i}")
            if is_not_assigned(group_idx):
                pattern = get_pattern(group_idx)
                
                #find group with max similarity
                min_sim = float('inf')
                min_idx = -1
                for j in range((i+1), self.N//current_block_size):
                    other_group_idx = j*current_block_size

                    if is_not_assigned(other_group_idx):
                        other_pattern = get_pattern(other_group_idx)
                        tmp_sim = similarity(pattern, other_pattern)
                        if tmp_sim < min_sim:
                            min_sim = tmp_sim
                            min_idx = other_group_idx
                
                #merge group with max similarity
                
                assign_group(group_idx,current_group_name)
                assign_group(min_idx, current_group_name)
                #print(f"merged block {i}, true {induced_row_order[group_idx]}, with block {max_idx}, true {induced_row_order[max_idx]}")
                #print(new_group)
                current_group_name += 1
        
        current_block_size *= 2;
        if current_block_size >= target_block_size:
            return new_group
        else:
            return self.hierchical_blocking(target_block_size, similarity, current_block_size, new_group)

                        
        
    def fixed_size_blocking(self, block_size, merge_crit):
        group_name = -1;
        group_array = np.ones(self.N)*(-1)
        
        for row_idx in range(self.N):
            if group_array[row_idx] == -1:
                group_name += 1;
                group_size = 1
                group_array[row_idx] = group_name
                pattern = self.pos[row_idx]
                size_reached = False
                for other_row_idx in range(row_idx+1,self.N):
                    if size_reached:
                        break
                    if group_array[other_row_idx] == -1:
                        other_pattern =  self.pos[other_row_idx]
                        candidates = len([x for x in group_array[other_row_idx + 1:] if x == -1])
                        merge = merge_crit(self.M, self.density, candidates, pattern, block_size - group_size, other_pattern)
                        if merge:
                            group_size += 1;
                            candidates = -1;
                            group_array[other_row_idx] = group_name;
                            pattern = merge_patterns(pattern, other_pattern)
                            if group_size == block_size:
                                size_reached = True
        return group_array
    
    def blocking(self,grouping, sim = None, verbose = False):
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
            
            size = 0
            while new_idx < len(induced_row_order) and grouping[induced_row_order[new_idx]] == current_group:
                
                

                old_idx = induced_row_order[new_idx]
                row  = self.pos[old_idx]
                new_idx +=1
                size += 1
                
                old_zeros = len(np.setdiff1d(current_pattern, row, assume_unique = True))                    
                if size != 1:
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
            print(f"TRUE NONZEROS: {self.tot_nz()} FAKE NONZEROS : {fake_nz}, with AVG. in-block DENSITY: {density}")
            print("PRINTING BLOCK DISTRIBUTION: size -- blocks with that size")
        
            for num in sorted(blocks_distribution):
                print(f"{num} --> {blocks_distribution[num]}")
            
        return density, n_of_block_rows, nz_blocks

    
    
    def blocking_show(self,grouping, filename = "blocking_structure_example.txt"):
        induced_row_order = sorted(range(len(grouping)), key=lambda k: grouping[k])
        
        current_group = 0;
        current_pattern = [];
           
        
        with open(filename, "w") as outfile:
            new_idx = 0
            n_of_block_rows = 0
            total_block_height = 0
            nz_blocks = 0
            fake_nz = 0
            blocks_distribution = {}
            while new_idx < len(induced_row_order):
                old_idx = induced_row_order[new_idx]
                
                fake_nz_here = 0
                size = 0
                while new_idx < len(induced_row_order) and grouping[induced_row_order[new_idx]] == current_group:
                    
                    

                    old_idx = induced_row_order[new_idx]
                    row  = self.pos[old_idx]
                    new_idx +=1
                    size += 1
                    
                    
                    old_zeros = len(np.setdiff1d(current_pattern, row, assume_unique = True))                    
                    
                    if size != 1:
                        fake_nz_here += old_zeros + len(np.setdiff1d(row, current_pattern,assume_unique = True))
                    
                    printrow = [0 for _ in range(old_zeros)] + [1 for _ in range(len(row))]
                    
                    print(printrow, file = outfile)
        
                    current_pattern = np.union1d(current_pattern,row);
                
                if len(current_pattern) != 0:
                    fake_nz += fake_nz_here
                    n_of_block_rows += 1
                    total_block_height += size*len(current_pattern)
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
            print(f"BLOCK ROWS: {n_of_block_rows} of AVG. SIZE: {total_block_height/nz_blocks}", file = outfile)
            print(f"TRUE NONZEROS: {self.tot_nz()} FAKE NONZEROS : {fake_nz}, with AVG. in-block DENSITY: {density}", file = outfile)
            print("PRINTING BLOCK DISTRIBUTION: size -- blocks with that size", file = outfile)
            
            for num in sorted(blocks_distribution):
                print(f"{num} --> {blocks_distribution[num]}", file = outfile)
            
            
        return density, total_block_height, nz_blocks
    
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


def weighted_sim(size1, p1, p2, tau, use_size = False, relative_val = True, cosine = False):
    unsim = general_sim(size1, p1, p2, use_size, relative_val, cosine);
    return unsim <= tau 
    
def general_sim(size1,p1,p2, use_size = False, relative_val = True, cosine = False):
    if (len(p1) == 0) and (len(p2) == 0):
        return 0
    
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
    return unsim

def merge_prob(M, density,  candidates, pattern, target_size, row, safety_mult = 2):
    if len(pattern) == 0 and len(row) == 0:
        return True
    if len(row) == 0 or len(pattern) == 0:
        return False
    prob = prob_of_better(M, pattern, row , density())
    best_for_pattern =  int(prob*candidates) < target_size*safety_mult
    #print(pattern, row, candidates, prob, merge)
    return best_for_pattern 

def double_merge_prob(M, density, candidates, pattern, target_size, row, safety_mult = 2):
    if len(pattern) == 0 and len(row) == 0:
        return True
    if len(row) == 0 or len(pattern) == 0:
        return False
    prob = prob_of_better(M, pattern, row , density())
    prob2 = prob_of_better(M, row, pattern, density())
    final_prob = min(prob,prob2)
    return  int(final_prob*candidates) < target_size*2



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
    
    
    


def Hamming(v1,v2):
    return len(np.setdiff1d(v1,v2)) + len(np.setdiff1d(v2,v1))


def Jaccard(v1,v2):
    return Hamming(v1,v2)/len(np.union1d(v1,v2))

def Special(v1,v2):
    if len(v1) == 0 and len(v2) == 0:
        return 0
    if len(v1) == 0 or len(v2) == 0:
        return 1
    num = Hamming(v1,v2)
    #den = len(v1) + len(v2) - num
    den = len(v1) + len(v2)
    if den == 0:
        return 1
    return num/den

def gHamming(v1, n1, v2, n2, mode = 0):
    if mode == 0:
        a = n2
        b = n1
    elif mode == 1:
        a = n1
        b = n2
        
    return a*len(np.setdiff1d(v1,v2)) + b*len(np.setdiff1d(v2,v1))


def gJaccard(v1, n1, v2, n2, mode = 0):
    if len(v1) == 0 and len(v2) == 0:
        return 0
    if len(v1) == 0 or len(v2) == 0:
        return 1
    num = gHamming(v1,n1,v2,n2, mode)
    den = n1*len(v1) + n2*len(v2) + gHamming(v1,n1,v2,n2)
    return num/den

def gSpecial(v1, n1, v2, n2, mode = 0):
    if len(v1) == 0 and len(v2) == 0:
        return 0
    if len(v1) == 0 or len(v2) == 0:
        return 1
    num = gHamming(v1,n1,v2,n2, 1)
    den = n1*len(v1) + n2*len(v2) # + gHamming(v1,n1,v2,n2,1)
    return num/den
    
    
jaccard = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = True)
jaccard_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = True)

hamming = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = False)
hamming_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = False)

cosine = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = False, relative_val = False, cosine = True)
cosine_special = lambda x,y,z,tau : weighted_sim(x,y,z,tau,use_size = True, relative_val = False, cosine = True)

similarities = {
    "jaccard" : lambda v1,n1,v2,n2: gJaccard(v1,1,v2,1),
    "jaccard_special" : gJaccard,
    "hamming": hamming,
    "hamming_special": hamming_special,
    "cosine" : cosine,
    "cosine_special": cosine_special,
    "special": lambda v1,n1,v2,n2: gSpecial(v1,n1,v2,1)
    }

fixed_size_criteria = {
        "fixsize_single_prob" : merge_prob,
        "fixsize_double_prob" : double_merge_prob
        }