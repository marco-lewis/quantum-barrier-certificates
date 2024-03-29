import itertools

import numpy as np
import sympy as sym
from scipy.optimize import NonlinearConstraint

z0 = sym.Symbol("z0", complex = True)
z1 = sym.Symbol("z1", complex = True)
z2 = sym.Symbol("z2", complex = True)
z3 = sym.Symbol("z3", complex = True)

def qubit_simple_constraints(zero_bounds, one_bounds):
    def f(x): return [x[0]**2 + x[1]**2,
                      x[2]**2 + x[3]**2,
                      x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2]
    return [NonlinearConstraint(f, [zero_bounds[0], one_bounds[0], 1], [zero_bounds[1], one_bounds[1], 1])]

close_to_zero = qubit_simple_constraints((0.9,1),(0,0.1))
close_to_one = qubit_simple_constraints((0,0.1),(0.9,1))

def generate_term_powers(k, num_complex_variables):
    x = range(k + 1)
    l = [list(p) for p in itertools.product(x, repeat=2*num_complex_variables)]
    less_than_k = filter(lambda c: sum(c) <= k, l)
    return list(less_than_k)

def round_sympy_expr(expr, prec=2):
    ret_expr = expr
    for a in sym.preorder_traversal(expr):
        if isinstance(a, sym.Float):
            ret_expr = ret_expr.subs(a, round(a, prec))
    return ret_expr

def get_real_from_model(model, var, prec=2):
    return round(float(model[var].numerator_as_long())/float(model[var].denominator_as_long()),prec)