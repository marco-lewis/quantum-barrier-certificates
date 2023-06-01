from __future__ import annotations
from copy import deepcopy
from typing import List

from src.utils import *

import numpy as np
import sympy as sym

class FuncTerm():
    def __init__(self, coeff, var):
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
        
    def __add__(self, other : FuncTerm):
        if other.var == self.var: return FuncTerm(self.coeff + other.coeff, self.var)
        return FuncSum([self, other])
    
    def __mul__(self, other : FuncTerm):
        var = [s + t for s, t in zip(self.var, other.var)]
        return FuncTerm(self.coeff * other.coeff, var)
    
    def __repr__(self):
        return "FuncTerm(" + repr([self.coeff, self.var]) + ")"
    
    def __str__(self):
        return str([str(self.coeff), str(self.var)])
    
class FuncSum():
    def __init__(self, fterms : List[FuncTerm]):
        self.fterms = fterms
    
    def conj(self):
        return FuncSum([t.conj() for t in self.fterms])
        
    def get_coeffs(self):
        return [t.coeff for t in self.fterms]
    
    def get_sym_sum(self, model, prec = 2):
        s = 0
        if isinstance(model, np.ndarray):
            for term, m in zip(self.fterms, model):
                v = term.get_var_sym()
                s += m * v
        else:
            for term in self.fterms:
                v = term.get_var_sym()
                c = term.coeff
                cr = get_real_from_model(model, c.r, prec)
                ci = get_real_from_model(model, c.i, prec)
                s += (cr + 1j*ci) * v
        return s
        
    def __add__(self, other : FuncSum):
        fterms = deepcopy(self.fterms)
        for t in other.fterms:
            modified = False
            for st in fterms:
                if t.var == st.var and not modified:
                    st.coeff = st.coeff + t.coeff
                    modified = True
            if not modified: fterms.append(t)
        return FuncSum(fterms)
    
    def __sub__(self, other : FuncSum):
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
    def __init__(self, fsums : List[FuncSum]):
        self.fsums = fsums
        
    def conj(self):
        return FuncVec([s.conj() for s in self.fsums])
        
    def __add__(self, other : FuncVec):
        return FuncVec([s0 + s1 for s0, s1 in zip(self.fsums, other.fsums)])
    
    # Dot product
    def __mul__(self, other : FuncVec):
        fsums = [s0 * s1 for s0, s1 in zip(self.fsums, other.fsums)]
        out = FuncSum([])
        for s in fsums:
            out += s
        return out
    
    def __repr__(self):
        return "FuncVec(" + repr(self.fsums) + ")"
    
    def __str__(self):
        return str(self.fsums)