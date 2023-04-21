from copy import copy

from src.FuncClasses import *
from src.utils import *

def sym2sum(sym_expr, num_vars=1):
    poly_len = num_vars * 2
    powers = [0] * poly_len
    terms = sym_expr.args

    fterms = []
    for term in terms:
        term_powers = copy(powers)
        vars = []
        for arg in term.args:
            exp = 1
            if isinstance(arg, sym.Pow): exp, arg = arg.exp, arg.base
            vars.append(arg)
            if isinstance(arg, sym.Symbol): term_powers[int(arg.name[1:])] += exp
            if isinstance(arg, sym.conjugate): term_powers[num_vars + int(arg.args[0].name[1:])] += exp
    
        if vars == []: c = float(term)
        else: c = sym.Poly(term, vars).coeffs()[0]
        fterms.append(FuncTerm(c, term_powers))
    
    return FuncSum(fterms)

def symlist2vec(sym_list, num_vars=1):
    return FuncVec([sym2sum(sym_expr, num_vars=num_vars) for sym_expr in sym_list])