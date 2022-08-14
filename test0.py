# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 10:47:51 2022

@author: Paolo
"""

import csr
import similarities as sims
import cmat_reorderings as blk
import experiments as exps
import numpy as np





print("meow")
cmat = csr.make_random_CSR(10,10,0.1);
taus = np.linspace(0, 1, 30);

block_size = 2
#results = exps.run_experiment(cmat, taus, block_size, opt_dict = HC)