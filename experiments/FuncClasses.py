import cmath
from copy import copy, deepcopy
import itertools

from complex import *
from ComplexVector import *
from FuncClasses import *
from utils import *

# import cvxpy as cp
import numpy as np
import sympy as sym
from z3 import *
# import scipy as sp

class Coeff():
    def __init__(self, nums, z3_vars):
        if len(nums) != len(z3_vars): raise Exception()
        self.nums = nums
        self.z3_var = z3_vars
    
    def __repr__(self):
        return "Coeff(" + _ + ")"

class FuncTerm():
    def __init__(self, coeff, var, coeff_factor=None):
        self.coeff = coeff
        self.var = var
        
    def conj(self):
        n = len(self.var)//2
        if isinstance(self.coeff, complex): return FuncTerm(np.conj(self.coeff), self.var[n:2*n] + self.var[:n])
        return FuncTerm(self.coeff.conj(), self.var[n:2*n] + self.var[:n])
    
    def get_var_sym(self):
        n = len(self.var)//2
        v = 1
        for i in range(n):
            z = sym.Symbol('z'+str(i), complex=True)
            v *= z**self.var[i] * np.conj(z)**self.var[i+n]
        return v
        
    def __add__(self, other):
        if other.var == self.var: return FuncTerm(self.coeff + other.coeff, self.var)
        return FuncSum([self, other])
    
    def __mul__(self, other):
        var = [s + t for s, t in zip(self.var, other.var)]
        return FuncTerm(self.coeff * other.coeff, var)
    
    def __repr__(self):
        return "FuncTerm(" + repr([self.coeff, self.var]) + ")"
    
    def __str__(self):
        return str([str(self.coeff), str(self.var)])
    
class FuncSum():
    def __init__(self, fterms):
        self.fterms = fterms
        
    def apply(self, nums):
        nums = nums + list(np.conj(nums))
        nums = [ComplexExpr(n.real, n.imag) for n in nums]
        s = 0
        for term in self.fterms:
            m = 1
            for n, v in zip(nums, term.var): m *= (n**v)
            s += (term.coeff * m)
        return s
    
    def conj(self):
        return FuncSum([t.conj() for t in self.fterms])
        
    def get_coeffs(self):
        return [t.coeff for t in self.fterms]
    
    def get_sym_sum(self, model, prec = 2):
        s = 0
        if isinstance(model, np.ndarray):
            for term, m in zip(self.fterms, model):
                v = term.get_var_sym()
                s += round(m, prec) * v
        else:
            for term in self.fterms:
                v = term.get_var_sym()
                c = term.coeff
                cr = get_real_from_model(model, c.r, prec)
                ci = get_real_from_model(model, c.i, prec)
                s += (cr + 1j*ci) * v
        return s
        
    def __add__(self, other):
        fterms = deepcopy(self.fterms)
        for t in other.fterms:
            modified = False
            for st in fterms:
                if t.var == st.var and not modified:
                    st.coeff = st.coeff + t.coeff
                    modified = True
            if not modified: fterms.append(t)
        return FuncSum(fterms)
    
    def __sub__(self, other):
        other = deepcopy(other)
        for t in other.fterms:
            t.coeff = 0- t.coeff
        return self + other
    
    def __mul__(self, other):
        return FuncSum([s*t for s, t in itertools.product(self.fterms, other.fterms)])
    
    def __repr__(self):
        return "FuncSum(" + repr(self.fterms) + ")"
    
    def __str__(self):
        return str(self.fterms)

class FuncVec():
    def __init__(self, fsums):
        self.fsums = fsums
        
    def conj(self):
        return FuncVec([s.conj() for s in self.fsums])
        
    def __add__(self, other):
        return FuncVec([s0 + s1 for s0, s1 in zip(self.fsums, other.fsums)])
    
    # Dot product
    def __mul__(self, other):
        fsums = [s0 * s1 for s0, s1 in zip(self.fsums, other.fsums)]
        out = FuncSum([])
        for s in fsums:
            out += s
        return out
    
    def __repr__(self):
        return "FuncVec(" + repr(self.fsums) + ")"
    
    def __str__(self):
        return str(self.fsums)
    
def diff(funcsum, var_loc):
    new_fsum = []
    for fterm in funcsum.fterms:
        t = deepcopy(fterm.var)
        if t[var_loc] != 0:
            c = fterm.coeff * t[var_loc]
            t[var_loc] -= 1
            new_fsum.append(FuncTerm(c, t))
    return FuncSum(new_fsum)