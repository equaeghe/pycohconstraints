#!/usr/bin/python

from timeit import timeit
from csv import writer
from pycohconstraints.support import randintlists
from pycohconstraints.coh import get_coh_hrep_via_vreps, get_coh_hrep_as_coeff_vrep

timeitonce = lambda stmt, setup: timeit(stmt, setup, number=1)

C1 = get_coh_hrep_via_vreps
C4 = get_coh_hrep_as_coeff_vrep

datafile = open('coh_cmp-128+.dat', 'a', 0)
datawriter = writer(datafile, delimiter='\t')

gamblenum = 5

fix_scope = "from __main__ import "

datawriter.writerow(['states', 'gambles', 'sparsty', 'fC1t', 'fC4t'])

for statenum in [2048 * (2 ** i) for i in xrange(4)]:
    for k in xrange(25):
        K, sparsity = randintlists(gamblenum, statenum)
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
        datawriter.writerow([statenum, gamblenum, sparsity, fC1t, fC4t])
