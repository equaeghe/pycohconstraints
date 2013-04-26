#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    """ Determines if r is a number ('int' or 'float'), returns a 'bool'. """
    return(isinstance(r, int) or isinstance(r, float))


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
	    else:
		print(n)
    except ValueError:
        exit()
