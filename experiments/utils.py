import itertools

import numpy as np
import sympy as sym
from scipy.linalg import expm
from scipy.optimize import NonlinearConstraint


isqrt2 = 1/np.sqrt(2)
had = np.array([[isqrt2, isqrt2], [isqrt2, -isqrt2]])

trange = np.arange(0, 2*np.pi, 0.01)

X = np.array([[0,1],[1,0]])
Y = np.array([[0,1j],[-1j,0]])
Z = np.array([[1,0],[0,-1]])

def print_barrier(barrier): print("Generated barrier: ", barrier, "\n")

def qubit_simple_constraints(zero_bounds, one_bounds):
    def f(x): return [x[0]**2 + x[1]**2,
                      x[2]**2 + x[3]**2,
                      x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2]
    return [NonlinearConstraint(f, [zero_bounds[0], one_bounds[0], 1], [zero_bounds[1], one_bounds[1], 1])]

close_to_zero = qubit_simple_constraints((0.9,1),(0,0.1))
close_to_one = qubit_simple_constraints((0,0.1),(0.9,1))

def ham(op):
    return lambda t: expm(-1j *op * t)

def had_ham(t):
    return expm(-1j * had * t)

def generate_term_powers(k, num_complex_variables):
    x = range(k + 1)
    l = [list(p) for p in itertools.product(x, repeat=2*num_complex_variables)]
    less_than_k = filter(lambda c: sum(c) <= k, l)
    return list(less_than_k)

def round_sympy_expr(expr):
    ret_expr = expr
    for a in sym.preorder_traversal(expr):
        if isinstance(a, sym.Float):
            ret_expr = ret_expr.subs(a, round(a, 2))
    return ret_expr