#!/usr/bin/python

from timeit import timeit
from csv import writer
from pycohconstraints.support import randintlists
from pycohconstraints.coh import get_coh_hrep_via_vreps, get_coh_hrep_as_coeff_vrep

timeitonce = lambda stmt, setup: timeit(stmt, setup, number=1)

C1 = get_coh_hrep_via_vreps
C4 = get_coh_hrep_as_coeff_vrep

datafile = open('coh_cmp.dat', 'w', 0)
datawriter = writer(datafile, delimiter='\t')

gamblenum = 5

fix_scope = "from __main__ import "

datawriter.writerow(['states', 'gambles', 'sparsty',
                     'rC1t', 'rC4t', 'fC1t', 'fC4t'])

for statenum in [4 * (2 ** i) for i in xrange(9)]:
    for k in xrange(100):
        K, sparsity = randintlists(gamblenum, statenum)
        try:
            rC1t = '{:.2g}'.format(timeitonce("C1(K)",
                                              fix_scope + "K, C1"))
        except:
            rC1t = 'X'
        try:
            rC4t = '{:.2g}'.format(timeitonce("C4(K)",
                                              fix_scope + "K, C4"))
        except:
            rC4t = 'X'
        try:
            fC1t = '{:.2g}'.format(timeitonce("C1(K, num_type='float')",
                                              fix_scope + "K, C1"))
        except:
            fC1t = 'X'
        try:
            fC4t = '{:.2g}'.format(timeitonce("C4(K, num_type='float')",
                                              fix_scope + "K, C4"))
        except:
            fC1t = 'X'
        datawriter.writerow([statenum, gamblenum, sparsity,
                             rC1t, rC4t, fC1t, fC4t])
