#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Generate the coherence constraints for a given set of gambles.

Calling constraints.py without an argument or with incorrect ones will
print its usage instructions to standard output.  These instructions can
be found towards the end of this file in the INSTRUCTIONS string
literal.

The script can also be used as a python module.  The only functions of
interest to the typical user are the gamble generating ones at the
beginning of this file.  The main program is at the very end.

"""

from __future__ import print_function
from itertools import chain

# Functions for generating gambles, sets of gambles, and also a set of
# values.
#
# These functions are both used internally and meant for a user who
# wants to generate a suitable pickled set of gambles to provide as
# input to the script.

TOL = 1e-10  # Choose wisely, it is used to distinguish a float from 0.

from itertools import product, repeat, permutations

def pspace(n):
    """ Generate the constant-1-valued gamble (length-n 'tuple'). """
    return(tuple(repeat(1, n)))

def singletons(n):
    """ Generate the 'set' of 1-once, 0-else gambles (length-n 'tuple'). """
    perms = permutations(list(repeat(0, n-1)) + [1], n)
    return(set(perms))

def singleton_complements(n):
    """ Generate the 'set' of 0-once, 1-else gambles (length-n 'tuple'). """
    perms = permutations(list(repeat(1, n-1)) + [0], n)
    return(set(perms))

def nonconstant_indicators(n):
    """ Generate the 'set' of nontrivial 0/1-valued gambles (length-n 'tuple'). """
    prods = product([0, 1], repeat=n)
    K = filter(lambda g: max(g) == 1 and min(g) == 0, prods)
    return(set(K))

def values_based_gambles(n, vals):
    """ Generate a 'set' of vals-valued gambles (length-n 'tuple').

    The set generated is the subset of the set of all vals-valued
    gambles that with a maximum of 1 and minimum of 0.  So to avoid
    generating the empty set, it should hold that {0,1} <= vals < [0,1].

    """
    L = product(vals, repeat=n)
    K = filter(lambda g: max(g) == 1 and min(g) == 0, L)
    return(set(K))

def generate_values(k):
    """ Generate the 'set' of k+1 values 0, 1/k, 2/k, ..., (k-1)/k, 1. """
    vals = (i/float(k) for i in range(0, k+1))
    return(set(vals))

def gamble_compare(f, g):
    """ Order two gambles inversely according to the sum of their values """
    diff = sum(f)-sum(g)
    if abs(diff) < TOL:
       return(0)
    if diff > 0:
        return(-1)
    else:
        return(1)

# Functions that generate the coherence constraints
#
# A constraint $\sumprod{\lambda}{\lp}{\gambles}\leq\mu$ is written
# as $0\leq\mu-\sumprod{-\lambda}{\lp}{\gambles}$ and as such
# encoded in a dictionary; we store all nonzero coefficients as
# illustrated by the following example:
# #{'mu': 1, g: -1}# corresponds to $\lp g\leq1$.
#
# There are two main functions, the first self-contained, short one
# generates the lower bound constraints of the coherence criterion,
# the second one generates all the other ones.  It is not
# self-contained: to make the things that go on digestible, it relies
# on a whole host of other functions (some non-self-contained
# themselves), each performing a specific function.  We will introduce
# further subdivisions, introduced by a one line comment.
#
# These functions are only meant for internal use.

def lower_bound_constraints(K):
    """ Generates an 'iterable' of lower bound constraints (a 'dict').

    The lower bound constraints are P(g) >= min(g) = 0.
    As this function is only meant for internal use,
    it can be assumed that the minimum is indeed 0.

    """
    L = [{'mu': 0, g: 1} for g in K]
    return(L)

def coherence_constraints(n, K):
    """ Generates an 'iterable' of coherence constraints (a 'dict').

    The constraints are generated as follows.  From K, we generate a
    list Ns of subsets N of gambles.  Then, in successive steps, we
    transform this N into a list of constraints.  The first step is to
    generate the matrix A of the linear system  A*lambda = b  that will
    be used to calculate the coefficients lambda.  Then we calculate
    the common support S of the gambles in N, which is used to filter
    out those N for which the number of gambles is too large with
    respect to the support (for us to need to consider them).  In the
    last (big) step, we generate all of the constraints corresponding
    to each N.

    These coherence constraints do not include any bound constraints,
    just those that can be considered as inter-gamble comparisons.

    """
    Ns = gamble_selections(n, K)
    X = map(to_NA, Ns)
    map(NA_to_NSA, X)
    X = filter(large_enough_support, X)
    Y = map(lambda x: NSA_to_constraints(K, x), X)
    return(chain.from_iterable(Y))

# coherence\_constraints helper functions

from itertools import combinations

from numpy import array, column_stack, sum, flatnonzero, mat
from numpy.linalg import qr, solve

def gamble_selections(n, K):
    """ Generate an 'iterable' of coherence constraint gamble 'list'.

    We take all combinations of 2 to n gambles.
    The bounds on the numbers of gambles follows from including the
    lower bound constraints and from linear independence.

    """
    N = chain.from_iterable(combinations(K, i) for i in range(2,n+1))
    return(N)

def to_NA(N):
    """ Generate a 'dict' with N and the matrix A built from it. """
    return({'N': N, 'A': column_stack(N)})

def NA_to_NSA(x):
    """ Add the support (vector of 'bool') of the gambles in N to x. """
    x['S'] = sum(x['A'], axis=1) != 0

def large_enough_support(x):
    """ Check if N's support S is large enough, returns 'bool'. """
    return(len(x['N']) <= len(flatnonzero(x['S'])))

def NSA_to_constraints(K, x):
    """ Generate the 'iterable' of constraints ( 'dict') corresponding to N.

    We first read out all the data in the dict.  Then we do a
    QR-factorization of A restricted to the support of N.  R is first
    used to see if the gambles in N are linearly independent (otherwise
    we do not need to consider this N).  Next we select the potential
    vectors b of the linear system  A*lambda = b, based on whether A
    and b have common support and gather them in a matrix B.  Then,
    using Q and R, we solve the multiple right-hand side linear system
    A*L = B.  After this, we reorganize the one x to multiple ones, one
    for each column lambda in L and b in B.  These are then checked to
    see if they indeed provide solutions to  A*lamba = b  (because we
    work with the thin QR-factorization, this is not always guaranteed)
    and satisfy the sign conditions.  The ones that pass the test get
    converted to constraints and are returned as function output in
    this form.

    """
    A = x['A']
    S = x['S']
    N = x['N']
    Q, R = qr(A[S])
    if lin_dep(R):
        x.clear()
        return([])
    coN = K  | set([pspace(n)]) - set(N)
    listB = list(filter(lambda g: common_support(S, g), coN))
    if len(listB) == 0:
        x.clear()
        return([])
    B = column_stack(listB)
    L = array(solve(R, mat(Q).T * mat(B[S])))  # array to allow indexation
    y = [x_to_xplus(x, B[:,i], L[:,i]) for i in range(0, len(listB))]
    sely = filter(lambda x: check_sgn(x) and check_sol(x), y)
    map(x_to_constraint, sely)
    return(sely)

# NSA\_to\_constraints helper functions

from numpy import amin, abs, diag, all

def common_support(S, g):
    """ Check if g has support S, returns 'bool'. """
    return(all((array(g) != 0) == S))

def lin_dep(R):
    """ Check if R (from a QR-factorization) has full rank, returns 'bool'. """
    r = amin(abs(diag(R)))
    return(r < TOL)

def x_to_xplus(x, b, lmbda):
    """ Return a copy of x ( 'dict') augmented with b and lmbda. """
    xplus = x.copy()
    xplus['b'] = tuple(b)
    xplus['lmbda'] = lmbda
    return(xplus)

def check_sol(x):
    """ Check if the QR-solution in x is a true solution, returns 'bool'. """
    b = mat(x['b'])
    A = mat(x['A'])
    lmbda = mat(x['lmbda'])
    return(all(abs(b.T - A*lmbda.T) < TOL))

def check_sgn(x):
    """ Checks if lmbda in x has correct signs, returns 'bool'.

    There are two cases to consider: if b is the constant-1-valued
    gamble, then lmbda may contain one negative component, otherwise
    all components must be strictly positive.  In neither case may a
    component of lmbda be zero.

    """
    b = x['b']
    lmbda = x['lmbda']
    zz = all(abs(lmbda) > TOL)  # zero zeroes
    non = len(flatnonzero(lmbda< TOL))  # number of nonpositives
    if set(b) == set([1]):
        return(zz and non <= 1)
    else:
        return(zz and non == 0)

def x_to_constraint(x):
    """ Modifies x into a constraint.

    We clean up objects in x that are not needed any more and add the
    necessary (gamble,coefficient)-pairs.

    """
    del x['S']
    del x['A']
    if all(x['b']) == 1:
        x['mu'] = 1
    else:
        x[x['b']] = 1
    del x['b']
    lmbda = x['lmbda']
    del x['lmbda']
    N = x['N']
    del x['N']
    for i in range(0,len(lmbda)):
        x[N[i]] = -lmbda[i]


# Functions used to validate user input.
#
# These functions are only meant for internal use.

def check_input(K):
    """ Checks if K is a valid 'set' of gambles, returns their length or 0.

    A set of gambles is valid if (i) it is non-empty and if
    (ii) it consists of valid gambles that have the same length.

    """
    if isinstance(K, set):
        if len(K) > 0:
            lensK = map(check_gamble, K)
            n = max(lensK)
            if min(lensK) == n:
                return(n)
    return(0)

def check_gamble(g):
    """ Checks if g is a valid gamble, returns g's length or 0.

    A gamble is valid if (i) it is a non-zero-length 'tuple' consisting
    of 'int' and/or 'float' and if (ii) it's maximum is 1 and minimum
    is 0.

    """
    if isinstance(g, tuple):
        if all(map(is_num, g)):
            if max(g) == 1 and min(g) == 0:
                return(len(g))
    return(0)

def is_num(r):
    """ Determines if r is a number ( 'int' or 'float'), returns a 'bool'. """
    return(isinstance(r, int) or isinstance(r, float))


# Functions used to generate the output of the script.
# The output is "Polyhedra H-Format", see the cddlib reference manual
# for more information (to be found on-line and distributed with the
# cdd program).
#
# These functions are only meant for internal use.

def print_constraints(n, K, X):
    """ Print the constraints, preceded by some useful information.

    As useful information, we provide the cardinality of the
    possibility space, the set of gambles used, ordered
    (which is important to interpret the constraints), and
    a brief explanation on how to interpret the constraints.
    Then, the list of constraints is given (between begin/end).

    """
    print("* Cardinality:", n)
    print("* The set of gambles, and their order:")
    print("*     K =", end=' ')
    print_gambles(K)
    print(EXPLANATION)
    print("H-representation")
    print("begin")
    print(len(X), len(K) + 1, "integer")
    for x in X:
      print_constraint(K, x)
    print("end")

def print_gambles(K):
    """ Print out the 'list' of gambles ( 'tuple').

    Before actually printing out the gambles, we convert their values
    (these are 'float'), to fractions to improve legibility.

    """
    fracK = [map(frac_ify, g) for g in K]
    print('{', end='')
    print_gamble(fracK.pop(0))
    for g in fracK:
        print(', ', end='')
        print_gamble(g)
    print('}')

def print_gamble(g):
    """ Print out a gamble ( 'tuple'). """
    print('(', g.pop(0), sep='', end='')
    for r in g:
        print(',', r, sep=' ', end='')
    print(')', end='')

EXPLANATION = """\
* A constraint is of the form
*     sum_{g in K} lambda_g P(g) <= mu.
* It is listed below as
*     mu -lambda_f -lambda_g -lambda_h ...\
"""

def print_constraint(K, x):
    """ Print a constraint.

    To print out the constraint x (a 'dict' of
    (gamble: coefficient)-pairs), of which only the non-zero
    coefficients have been kept up until now, we first embed it in a
    full constraint, using an all-zero template created using their
    'list' of gambles K.

    """
    e = dict.fromkeys(K, 0)
    e['mu'] = 0
    e.update(x)
    print(e['mu'], end=' ')
    for g in K:
      print(e[g], end=' ')
    print() # newline

INSTRUCTIONS = """
 constraints.py generates the coherence constraints
 for a given set of gambles in Polyhedra H-format
 (for details about this format, see the cddlib reference
  manual, available on the internet).

 This set of gambles has to be saved as a pickled set
 of tuples of uniform (but arbirary) length containing
 only (essentially) rational numbers between and
 including zero and one.

 An example will make it clear how to do this.
 First, enter the python interpreter:

  $ python
  >>> from pickle import *
  >>> output = open('K.pkl', 'wb')
  >>> mysetofgambles = set([(0,1,1),(1,.5,0),(1,0,1/3.)])
  >>> mysetofgambles
  set([(1, 0.5, 0), (0, 1, 1), (1, 0, 0.33333333333331)])
  >>> dump(mysetofgambles, output)
  >>> output.close()

 Typing Ctrl-D exits the python interpreter.
 constraints.py also acts as a python module containing
 some useful gamble-generating functions:

  >>> from constraints import *
  >>> nonconstant_indicators(2)
  set([(0, 1), (1, 0)])
  >>> singletons(3)
  set([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
  >>> singleton_complements(3)
  set([(0, 1, 1), (1, 1, 0), (1, 0, 1)])

 Peruse the first part of constraints.py for
 an overview of these useful functions.

 Now use the pickled set we created:

  $ ./constraints.py K.pkl
* Cardinality: 3
* The set of gambles, and their order:
*     K = {(0, 1, 1), (1, 1/2, 0), (1, 0, 1/3), \
(1, 0, 0), (0, 1, 0), (0, 0, 1)}
* A constraint is of the form
*     sum_{g in K} lambda_g P(g) <= mu.
* It is listed below as
*     mu -lambda_f -lambda_g -lambda_h ...
H-representation
begin
23 7 integer
0 0 1 0 0 0 0
0 1 0 0 0 0 0
0 0 0 0 1 0 0
0 0 0 0 0 1 0
0 0 0 1 0 0 0
0 0 0 0 0 0 1
1 -1 0 0 -1 0 0
0 0 2 0 -2 -1 0
0 0 0 3 -3 0 -1
0 1 0 0 0 -1 -1
2 -2 -2 0 0 1 0
5 -4 -2 -3 0 0 0
2 -1 -2 0 0 0 -1
1 0 -2 -3 4 0 0
1 0 -2 0 1 0 -1
1 0 2 -3 0 -2 0
2 0 -2 0 0 -1 -2
3 0 -6 3 0 0 -4
3 -2 0 -3 0 -1 0
3 -3 0 -3 0 0 1
1 0 0 -3 2 -1 0
1 0 0 0 -1 -1 -1
3 0 0 -3 0 -3 -2
end


 Note that the set of gambles used is larger than the
 pickled one we provided.  The set of singleton
 indicators is always added by the script.

 This output can be directly used as input to
 vertex-enumeration programs such as cdd and lrs.
"""


# Functions to transform constraints expressed using 'float' to ones
# expressed using 'int'.  Includes a function to transform 'float'
# to 'Fraction'.
#
# These functions are only meant for internal use.

from fractions import Fraction, gcd

MAXDENOM = pow(10,10)   # Choose wisely; to convert 'float' to 'Fraction'

def int_ify(x):
    """ Transform a 'float'-based constraint to an 'int'-based one.

    This function first converts all the 'float' to 'Fraction',
    then computes the least common multiple of the appearing
    denominators, by which all 'Fraction' are then multiplied.

    """
    for key in x:
        x[key] = frac_ify(x[key])
    denomlcm = reduce(lcm, [r.denominator for r in x.itervalues()])
    for key in x:
        x[key] = x[key].numerator * (denomlcm / x[key].denominator)

def frac_ify(r):
    """ Convert a 'float' to a 'Fraction' (max. denominator is MAXDENOM). """
    q = Fraction.from_float(r)
    return(q.limit_denominator(MAXDENOM))

def lcm(a,b):
    """ Calculate the least common multiple (an 'int') of two 'int'. """
    return(abs(a*b)/gcd(a,b))


# The 'main' part of the script.
#
# It first reads and validates the script's argument, which should be
# the name of a file containing a pickled 'set' of gambles.  These are
# equal-length 'tuple' with 'int' and/or 'float' components that have
# a maximum of 1 and minumum of 0 and are essentially rationals with a
# denominator no greater than MAXDENOM.  The intent is for the user to
# generate and pickle this object in the python interpreter.
#
# If the input does not validate, we show the instructions and exit.
#
# If the input validates, we add some essential gambles to the
# user-supplied set, generate the coherence constraints---lower bound
# constraints separately---, recast the constraints as integer-valued
# ones, and print them out in a format suitable for vertex enumerating
# programs such as cdd and lrs.

if __name__ == "__main__":

    from sys import argv, exit
    from pickle import load

    try:
        if len(argv) != 2:
            raise ValueError
        else:
            userKpkl = open(argv[1], 'rb')
            userK = load(userKpkl)
            userKpkl.close()
            n = check_input(userK)
            if n == 0:
                raise ValueError
    except ValueError:
        print(INSTRUCTIONS)
        exit()

    # Singletons are essential for the finite coherence constraint
    # generation algorithm.
    K = singletons(n) | userK

    B = lower_bound_constraints(K)
    C = coherence_constraints(n, K)
    X = list(chain(B, C))
    map(int_ify, X)

    print_constraints(n, sorted(list(K), gamble_compare), X)  # order matters
