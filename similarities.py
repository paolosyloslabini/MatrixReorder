# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 13:30:23 2022

@author: Paolo
"""

import numpy as np
import math


def get_block_pattern(p, block_size):
    res = []
    for elem in p:
        new = math.floor(elem/block_size)
        if (not res) or res[-1] != new:
            res.append(new)
    return res

def merge_patterns(p1,p2):
    i = 0
    j = 0
    res = []
    while (i < len(p1) and j < len(p2)):
        if p1[i] < p2[j]:
            res.append(p1[i])
            i += 1;
        elif p1[i] > p2[j]:
            res.append(p2[j])
            j += 1;
        else:
            res.append(p1[i]);
            i += 1;
            j += 1;
    
    if i < len(p1):
        res.extend(p1[i:])
    if j < len(p2):
        res.extend(p2[j:])
        
    return res
    
def HammingGeneral(p1,p2,block_size = 1, s1 = 1,s2 = 1, count_zeros = True):
  if (len(p1) == 0 and len(p2) == 0): return 0;

  i = 0;
  j = 0;
  count_1 = 0;
  count_2 = 0;

  blockpos_1 = lambda x: math.floor(p1[x]/block_size) 
  blockpos_2 = lambda x: math.floor(p2[x]/block_size)

  while (i < len(p1) and j < len(p2)):
    
    pos_A = blockpos_1(i);
    pos_B = blockpos_2(j);
    
    if (pos_A < pos_B):
      count_1 += 1
      i += 1
      while(i < len(p1) and blockpos_1(i) == pos_A): 
          i += 1;
      
    elif (pos_A > pos_B):
      count_2 += 1
      j += 1
      while(j < len(p2) and blockpos_2(j) == pos_B):
          j += 1;
          
    else:
      while(i < len(p1) and blockpos_1(i) == pos_A):
          i += 1;
      while(j < len(p2) and blockpos_2(j) == pos_B):
          j += 1;
          
  while (i < len(p1)):
    pos_A = blockpos_1(i);
    i += 1;
    count_1 += 1;
    while(i < len(p1) and blockpos_1(i) == pos_A):
        i += 1;
        
  while (j < len(p2)):
    pos_B = blockpos_2(j);
    j += 1;
    count_1 += 1;
    while(j < len(p2) and blockpos_2(j) == pos_B):
        j += 1;

  if count_zeros: 
      return count_1*s2 + count_2*s1
  else:
      return count_1*s1 + count_2*s2

    
def JaccardGeneral(p1,p2,block_size = 1, s1 = 1,s2 = 1, count_zeros = True):
    if (len(p1) == 0 and len(p2) == 0): return 0;
        
    h = HammingGeneral(p1,p2,block_size, s1,s2, count_zeros)
    return 2*h/(len(p1)*s1 + len(p2)*s2 + h)



Jaccard = lambda p1,p2,block_size : JaccardGeneral(p1,p2,block_size)
Hamming = lambda p1,p2,block_size : HammingGeneral(p1,p2,block_size)


