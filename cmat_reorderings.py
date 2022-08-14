# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 13:26:44 2022

@author: Paolo
"""
import numpy as np
import csr
from similarities import *
    
    

def IterativeClusteringSimple(cmat, tau, dist_func):
    
    group_name = -1;
    group_array = np.ones(cmat.N)*(-1)
    for row_idx in range(cmat.N):
        if group_array[row_idx] == -1:
            group_name += 1;
            group_size = 1
            group_array[row_idx] = group_name
            
            for other_row_idx in range(row_idx+1,cmat.N):
                if group_array[other_row_idx] == -1:
                    if dist_func(cmat.pos[row_idx], cmat.pos[other_row_idx],1,1) <= tau:
                        group_size += 1;
                        group_array[other_row_idx] = group_name;

    return group_array


def IterativeClusteringPattern(cmat, tau, dist_func):
    
    group_name = -1;
    group_array = np.ones(cmat.N)*(-1)
    for i in range(cmat.N):
        if group_array[i] == -1:
            group_name += 1;
            group_size = 1
            group_array[i] = group_name
            pattern = cmat.pos[i]
            
            for j in range(i+1,cmat.N):
                if group_array[j] == -1:
                    other_pattern =  cmat.pos[j]
                    if dist_func(pattern, other_pattern, group_size,1) <= tau:
                        group_size += 1;
                        group_array[j] = group_name;
                        pattern = merge_patterns(pattern, other_pattern)

    return group_array


def HierarchicalClustering(cmat, tau, dist_func):
    #dist_func(row_1, size_1, row_2, size_2)
    class Group():
         def __init__(self, cmat, row_idx):
             self.rows = [row_idx,]
             self.ID = row_idx
             self.cmat = cmat
             self.pattern = cmat.pos[row_idx]
             self.closest = -1
             self.closest_val = float('inf')
         
         def add(self,group):
             self.rows = group.rows + self.rows
             self.pattern = merge_patterns(self.pattern, group.pattern)


         def __eq__(self, other):
             if isinstance(other, Group):
                 return self.ID == other.ID
             return False
         
         def update_closest(self, sender):
             new_val = self.distance(sender)
             if new_val <= self.closest_val:
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
             d = dist_func(p1,p2,s1,s2)
             return d
             
         def find_closest(self, group_list):
             min_dist = float('inf')
             best_group = None
             for g in group_list:
                 if g.ID != self.ID:
                     tmp_dist = self.distance(g)
                     if tmp_dist <= min_dist:
                         min_dist = tmp_dist
                         best_group = g
             self.closest = best_group
             self.closest_val = min_dist
             if best_group == None:
                 print(f"DEBUG: grouplist size {len(group_list)}, group {self.ID} with pattern {self.pattern}, closest val {min_sim}")
   
    groups = [Group(cmat, row) for row in range(cmat.N)]
    
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
                
        
    grouping = np.ones(cmat.N)
    current = 0
    for g in groups:
        for row in g.rows:
            grouping[row] = current
        current += 1
    return grouping


blockings = {
    "IterativeClusteringSimple": IterativeClusteringSimple,
    "IterativeClusteringPattern" : IterativeClusteringPattern,
    "HierarchicalClustering": HierarchicalClustering
    }
