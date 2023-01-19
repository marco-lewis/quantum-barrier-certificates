from copy import deepcopy
import warnings

from FuncClasses import *
from utils import *

import numpy as np
import sympy as sym
from scipy.optimize import Bounds, minimize, linprog, OptimizeWarning

def diff_fsum(fsum, i):
    new_fsum = []
    for fterm in fsum.fterms:
        t = deepcopy(fterm.var)
        if t[i] != 0:
            c = fterm.coeff * t[i]
            t[i] -= 1
            new_fsum.append(FuncTerm(c, t))
    return FuncSum(new_fsum)

def scipy_find_b(barrier, states, f_vec, term_powers, prec=2, linprog_obj=[]):
    # Make appropriate conditions using representation
    dbdz = FuncVec([diff_fsum(barrier, i) for i in range(states)])
    dbdzconj = FuncVec([diff(barrier, i) for i in range(states, 2*states)])
    dbdt = (dbdz * f_vec) + (dbdzconj * f_vec.conj())
    barr_minus_conj = barrier - barrier.conj()
    
    # Create linprog matrices and vectors
    A_dbdt = [list(term.coeff.imag) for term in dbdt.fterms]
    b_dbdt = [0]*len(A_dbdt)
    
    A_conj = [list(term.coeff) for term in barr_minus_conj.fterms]
    b_conj = [0]*len(A_conj)
    
    A = np.array(A_dbdt + A_conj)
    b = np.array(b_dbdt + b_conj)
    # Minimize vector (min c^T x subject to ...)
    c = np.array([0]*len(term_powers))
    for i in linprog_obj: c[i] = 1
    
    # Negate A doesn't matter unless used for A_ub
    res = linprog(c,
                  A_eq=A,
                  b_eq=b,
                  bounds=(-1,1),
                 )
    return barrier.get_sym_sum(res.x, prec)

def scipy_find_constant(barrier_sym, states, init=[], prec=2):
    # Create objective function from Sympy barrier
    def obj(x):
        z = [x[2*i] + 1j*x[2*i + 1] for i in range(states)]
        b = barrier_sym
        for var_sym in barrier_sym.free_symbols:
            sym_num = int(str(var_sym)[1:])
            b = b.subs({var_sym: z[sym_num]})
        return -float(sym.re(b))

    # Bounds and guesses
    bounds = Bounds([-1]*2*states, [1]*2*states)
    guess = [0]*2*states
    
    res = minimize(obj,
                   guess,
                   method='trust-constr',
                   constraints=init,
                   options={'verbose': 0},
                   bounds=bounds,
                   hess=lambda x: np.zeros((2*states,))
                  )
    minimum = round(res.fun, prec)
    return minimum

def scipy_check_constant(c, barrier_sym, states, unsafe=[], prec=2):
    def obj(x):
        z = [x[2*i] + 1j*x[2*i + 1] for i in range(states)]
        b = barrier_sym
        for var_sym in barrier_sym.free_symbols:
            sym_num = int(str(var_sym)[1:])
            b = b.subs({var_sym: z[sym_num]})
        return float(sym.re(b))

    # Bounds and guesses
    bounds = Bounds([-1]*2*states, [1]*2*states)
    guess = [0]*2*states
    
    res = minimize(obj,
                   guess,
                   method='trust-constr',
                   constraints=unsafe,
                   options={'verbose': 0},
                   bounds=bounds,
                   hess=lambda x: np.zeros((2*states,))
                  )
    minimum = round(res.fun, prec)
    if -minimum >= c : raise Exception(str(barrier_sym) + "\nError: proposed barrier has part of unsafe in same contour as initial region")

# Setup coefficients as ndarray
def scipy_find_k_barrier(k, H, init=[], unsafe=[], linprog_obj = [], prec=2, verbose=False):
    z = -1j
    n = round(len(H))
    term_powers = generate_term_powers(k, n)
    coeff_num = len(term_powers)

    if verbose: print("Converting dynamical system...")
    sums = []
    for i in range(n):
        terms = []
        for j in range(n):
            t = [0]*(2*n)
            t[j] = 1
            t = tuple(t)
            terms.append(FuncTerm(z * H[i][j], t))
        sums.append(FuncSum(terms))
    f_vec = FuncVec(sums)
    if verbose: print("Dynamical system converted.")
    
    id_coeff = np.identity(coeff_num)
    barrier = FuncSum(list([FuncTerm(i, t) for i, t in zip(id_coeff, term_powers)]))
    if verbose: print("Finding polynomial...")
    warnings.simplefilter("ignore", OptimizeWarning)
    b = scipy_find_b(barrier, n, f_vec, term_powers, prec, linprog_obj)
    if verbose: print("Polynomial found: ", b)
    
    if verbose: print("Finding constant...")
    c = scipy_find_constant(b, n, init=init)
    if verbose: print("Checking...")
    scipy_check_constant(c, b, n, unsafe=unsafe)
    if verbose: print("Constant found: ", c)
    
    return round_sympy_expr(c + b)