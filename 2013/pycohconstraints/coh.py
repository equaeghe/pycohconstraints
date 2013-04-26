from pycohconstraints.support import transpose, diagonallists
from pycohconstraints.support import make_vrep, make_hrep
from pycohconstraints.asl import generate_asl_vrep, get_asl_lhs_rhs_as_coeff_vrep
from cdd import Matrix, RepType, Polyhedron
from itertools import chain

def generate_Sasl_pos_dirs(k, i):
    """Generate diagonal list of size `k` with 1s, apart from -1 at `i`"""
    pre = i
    post = k - i - 1
    return diagonallists([1] * pre + [-1] + [1] * post)

def generate_Sasl_neg_dirs(k, i):
    """Generate diagonal list of size `k` with -1s, apart from 1 at `i`"""
    pre = i
    post = k - i - 1
    return diagonallists([-1] * pre + [1] + [-1] * post)

def generate_Sasl_vrep(K, i, num_type='fraction'):
    """Generate the V-representation of avoiding S-sure loss

    See Equation (9) in my ISIPTA '13 paper “Characterizing coherence,
    correcting incoherence”. Here, S has negative component at index `i`.

    """
    deg_prevs = transpose(K)
    neg_dirs = generate_Sasl_neg_dirs(len(K), i)
    return make_vrep(deg_prevs, neg_dirs, num_type)

def get_coh_hrep_via_vreps(K, num_type='fraction'):
    """Compute a minimal H-representation for coherence

    See Procedure C.1 in my ISIPTA '13 paper “Characterizing coherence,
    correcting incoherence”.

    """
    vrep = generate_asl_vrep(K, num_type)
    hreplist = [Polyhedron(vrep).get_inequalities()]
    for i in xrange(len(K)):
        vrep = generate_Sasl_vrep(K, i, num_type)
        hreplist.append(Polyhedron(vrep).get_inequalities())
    constraints = ([list(t) for t in hrep[:]] for hrep in hreplist)
    constraintslist = list(chain.from_iterable(constraints))
    hrep = Matrix(constraintslist, number_type=num_type)
    hrep.rep_type = RepType.INEQUALITY
    hrep.canonicalize()
    return hrep

def get_coh_hrep_as_coeff_vrep(K, num_type='fraction'):
    """Compute a minimal H-representation for coherence

    See Procedure C.4 in my ISIPTA '13 paper “Characterizing coherence,
    correcting incoherence”.

    """
    k = len(K)
    lhs, rhs = get_asl_lhs_rhs_as_coeff_vrep(K)
    for i in xrange(k):
        S = generate_Sasl_pos_dirs(k, i)
        minusS = generate_Sasl_neg_dirs(k, i)
        coeff_one_lhs = transpose(K) + minusS
        coeff_one_rhs = ([1] * (len(coeff_one_lhs) - k)) + ([0] * k)
        coeff_one_hrep = make_hrep(coeff_one_lhs, coeff_one_rhs, num_type)
        coeff_one_vrep = Polyhedron(coeff_one_hrep).get_generators()
        lhs.extend(list(t[1:]) for t in coeff_one_vrep[:])
        rhs.extend(t[0] for t in coeff_one_vrep[:])
        coeff_zero_lhs = transpose(K) + minusS + S
        coeff_zero_rhs = ([0] * (len(coeff_zero_lhs) - k)) + ([1] * k)
        coeff_zero_hrep = make_hrep(coeff_zero_lhs, coeff_zero_rhs, num_type)
        coeff_zero_vrep = Polyhedron(coeff_zero_hrep).get_generators()
        lhs.extend(list(t[1:]) for t in coeff_zero_vrep[:])
        rhs.extend([0] * len(coeff_zero_vrep))
    pos_lhs = diagonallists([-1] * k)
    pos_rhs = ([0] * k)
    hrep = make_hrep(lhs + pos_lhs, rhs + pos_rhs, num_type)
    hrep.rep_type = RepType.INEQUALITY
    hrep.canonicalize()
    return hrep
