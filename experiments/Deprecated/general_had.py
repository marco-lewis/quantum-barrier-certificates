from typing import List
from utils_z3 import *

import numpy as np

def conj(z): return np.conj(z) if not(isinstance(z, ComplexExpr)) else z.conj()
    
H_hat = np.array([[1,1],[1,-1]])
H_hat_n = lambda n: np.kron(H_hat, H_hat_n(n-1)) if n > 1 else H_hat if n == 1 else 0

n = 2
N = 2**n
norm_sqrd = lambda z: z.r**2 + z.i**2 if isinstance(z, ComplexExpr) else z.real**2 + z.imag**2
sum_norm_sqrs = lambda z: sum([norm_sqrd(zi) for zi in z])

def evo(H :np.ndarray, z):
    dpsidt = [[0 -I * e for e in row] for row in H] if isinstance(z[0], ComplexExpr) else -1j*H
    n = len(H)
    return [sum([row[j] * z[j] for j in range(n)]) for row in dpsidt]

# Barrier definitions
z_const = 2
sum_term = lambda z: sum([z[0] * conj(zi) for zi in z[1:]]) + sum([conj(z[0]) * zi for zi in z[1:]])
barrier = lambda z: 1.16 - (z_const * z[0] * conj(z[0]) + sum_term(z))

def dbdz(z):
    diff = [-z_const * conj(z[0]) - np.sum([conj(zi) for zi in z[1:]])]
    for j in range(1, len(z)): diff.append(0 - conj(z[0]))
    return diff

def dbdzconj(z):
    diff_conj = [-z_const * z[0] - np.sum(z[1:])]
    for j in range(1, len(z)): diff_conj.append(0 - z[0])
    return diff_conj

# Defining barrier and Hamiltonian for n
b = barrier
ham = H_hat_n(n)
init = lambda z: Implies(norm_sqrd(z[0]) >= 0.9, And(b(z).r <= 0, b(z).i == 0))
# unsafe = lambda z: Implies(Or([norm_sqrd(zi) >= 0.9 for zi in z[1:]]), And(b(z).r > 0, b(z).i == 0))
unsafe = lambda z: Implies(norm_sqrd(z[1]) >= 0.9, And(b(z).r > 0, b(z).i == 0))
dbdt = lambda z: np.dot(dbdz(z), evo(ham, z)) + np.dot(dbdzconj(z), [conj(e) for e in evo(ham, z)])
evo_cond = lambda z: And(dbdt(z).r <= 0, dbdt(z).i == 0)
conds = [evo_cond, init, unsafe]
strs = ["Evolution", "Init", "Unsafe"]

if __name__ == "__main__":
    s = z3.Solver()
    z = ComplexVector('z', 2**n)
    s.add(sum_norm_sqrs(z) == 1)
    s.push()
    for cond, cond_str in zip(conds, strs):
        print(cond_str)
        s.add(Not(cond(z)))
        print("Checking...")
        sat = s.check()
        
        if sat == z3.sat:
            print("Found CEx")
            mz = [0] * 2**n
            m = s.model()
            for var in m:
                index, r_or_i = str(var)[3:].split('.')
                index = int(index)
                value = m[var()]
                if is_rational_value(value): v = value.as_fraction()
                elif is_algebraic_value(value):
                    r = value.approx(20)
                    v = float(r.numerator_as_long())/float(r.denominator_as_long())
                else: print(type(value))
                mz[index] +=  v if r_or_i == 'r' else v*1j
            print("Quantum State: ", mz)
            print("Norm: ",sum_norm_sqrs(mz))
            print("Barrier value: ", b(mz))
            print("Differential value: ", dbdt(mz))
        else: print(sat)

        print()
        s.pop()