from copy import deepcopy
import warnings

from src.FuncClasses import *
from src.utils import *

from scipy.optimize import Bounds, minimize, linprog, OptimizeWarning

def diff_fsum(funcsum, var_loc):
    new_fsum = []
    for fterm in funcsum.fterms:
        t = deepcopy(fterm.var)
        if t[var_loc] != 0:
            c = fterm.coeff * t[var_loc]
            t[var_loc] -= 1
            new_fsum.append(FuncTerm(c, t))
    return FuncSum(new_fsum)

# Step 1
def scipy_find_b(barrier : FuncSum, states, f_vec : FuncVec, term_powers, prec=2, obj_func_idxs=[]):
    # Make appropriate conditions using representation
    dbdz = FuncVec([diff_fsum(barrier, i) for i in range(states)])
    dbdzconj = FuncVec([diff_fsum(barrier, i) for i in range(states, 2*states)])
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
    if len(obj_func_idxs) == 0:
        c = np.ones(len(term_powers))
    else:
        c = np.zeros(len(term_powers))
        for idx in obj_func_idxs: c[idx] = 1

    # Negate A doesn't matter unless used for A_ub
    res = linprog(c,
                  A_eq=A,
                  b_eq=b,
                  bounds=(-1,1),
                 )
    return barrier.get_sym_sum(res.x, prec)

# Step 2
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

# Step 3
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

# Algorithm 1
def scipy_find_k_barrier(k, H, init=[], unsafe=[], prec=2, verbose=False, obj_func_idxs=[]):
    # Step 0 - setup data structures
    z = -1j
    n = round(len(H))
    term_powers = generate_term_powers(k, n)
    coeff_num = len(term_powers)

    if verbose:
        print("Step 0: Setting up data structures")
        print("Converting dynamical system...")
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
    
    if verbose: print("Creating barrier representation (B(a,z))...")
    id_coeff = np.identity(coeff_num)
    barrier = FuncSum(list([FuncTerm(i, t) for i, t in zip(id_coeff, term_powers)]))
    if verbose: print("Barrier created.")

    if verbose: print("Step 1: Finding polynomial (b)...")
    warnings.simplefilter("ignore", OptimizeWarning)
    b = scipy_find_b(barrier, n, f_vec, term_powers, prec, obj_func_idxs)
    if verbose: print("Polynomial found: b = " + str(b))
    
    if verbose: print("Step 2: Finding constant (c)...")
    c = scipy_find_constant(b, n, init=init, prec=prec)

    if verbose: print("Step 3: Checking constant (c)...")
    scipy_check_constant(c, b, n, unsafe=unsafe, prec=prec)
    if verbose: print("Constant found: c = " + str(c))
    barrier_out = round_sympy_expr(c + b)
    if verbose: print("Generated barrier: B(z) = b + c = " + str(barrier_out) + "\n")
    return barrier_out