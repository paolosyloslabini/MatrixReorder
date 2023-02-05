# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 17:22:21 2022

@author: Paolo
"""

import scipy.sparse as sprs
import matplotlib.pyplot as plt 
import numpy as np
Matrix=sprs.rand(10,10, density=0.2, format='csc')
Matrix[Matrix > 0] = 1


fig, ax = plt.subplots(1,1)

sparsity = 0.1


#make small lines

n = 10
m = 10
d = 0.2


x = np.random.rand(n,m)
x[x > d] = 0


plt.hlines(np.linspace(0.5, n + 0.5, n + 1),-1,m, color = "gray")
plt.vlines(np.linspace(0.5, m + 0.5, m + 1),-1,n, color = "gray")



hpart = np.array((-1,2,4,7,8.95))
vpart = np.array((-1,2,6,8.95))

plt.hlines(hpart + 0.5 ,-1,m, lw = 3, color = "red")
plt.vlines(vpart + 0.5 ,-1,n, lw = 3, color = "red")


ax.spy(x)

ax.axis('off')

plt.show()

