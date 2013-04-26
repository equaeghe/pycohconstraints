#!/usr/bin/python

from timeit import timeit
from csv import writer
from pycohconstraints.support import randintlists
from pycohconstraints.coh import get_coh_hrep_via_vreps

timeitonce = lambda stmt, setup: timeit(stmt, setup, number=1)

C1 = get_coh_hrep_via_vreps

datafile = open('coh_6states.dat', 'a', 0)
datawriter = writer(datafile, delimiter='\t')

statenum = 6

fix_scope = "from __main__ import "

datawriter.writerow(['states', 'gambles', 'sparsty', 'fC1t', 'fC1t'])

for gamblenum in [12]:
    for k in xrange(10):
        K, sparsity = randintlists(gamblenum, statenum)
        try:
            fC1t = '{:.2g}'.format(timeitonce("C1(K, num_type='float')",
                                              fix_scope + "K, C1"))
            rC1t = ''
        except:
            try:
                fC1t = 'X'
                rC1t = '{:.2g}'.format(timeitonce("C1(K)",
                                                  fix_scope + "K, C1"))
            except:
                rC1t = 'X'
        datawriter.writerow([statenum, gamblenum, sparsity, fC1t, rC1t])
