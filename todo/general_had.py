from typing import List
from experiments.todo.utils_z3 import *

import numpy as np

def conj(z): return np.conj(z) if not(isinstance(z, ComplexExpr)) else z.conj()
    
H_hat = np.array([[1,1],[1,-1]])
H_hat_n = lambda n: np.kron(H_hat, H_hat_n(n-1)) if n > 1 else H_hat if n == 1 else 0

n = 3
N = 2**n
norm_sqrd = lambda z: z.r**2 + z.i**2 if isinstance(z, ComplexExpr) else z.real**2 + z.imag**2
sum_norm_sqrs = lambda z: sum([norm_sqrd(zi) for zi in z])

def evo(H :np.ndarray, z):
    dpsidt = [[0 -I * e for e in row] for row in H] if isinstance(z[0], ComplexExpr) else -1j*H
    n = len(H)
    return [sum([row[j] * z[j] for j in range(n)]) for row in dpsidt]

# Barrier definitions
c = -1
barrier = lambda z: c + 1*z[0]**2 - (2/7)*z[0]*z[1] - (2/7)*z[0]*z[2] - (2/7)*z[0]*z[3] - (2/7)*z[0]*z[4] - (2/7)*z[0]*z[5] - (2/7)*z[0]*z[6] - (2/7)*z[0]*z[7] - 1*z[0]*conj(z[0]) - (1/9)*z[0]*conj(z[1]) - (1/9)*z[0]*conj(z[2]) - (1/9)*z[0]*conj(z[3]) - (1/9)*z[0]*conj(z[4]) - (1/9)*z[0]*conj(z[5]) - (1/9)*z[0]*conj(z[6]) - (1/9)*z[0]*conj(z[7]) - (1/7)*z[1]**2 - (2/7)*z[1]*z[2] - (2/7)*z[1]*z[3] - (2/7)*z[1]*z[4] - (2/7)*z[1]*z[5] - (2/7)*z[1]*z[6] - (2/7)*z[1]*z[7] - (1/9)*z[1]*conj(z[0]) - (1/9)*z[1]*conj(z[1]) - (1/9)*z[1]*conj(z[2]) - (1/9)*z[1]*conj(z[3]) - (1/9)*z[1]*conj(z[4]) - (1/9)*z[1]*conj(z[5]) - (1/9)*z[1]*conj(z[6]) - (1/9)*z[1]*conj(z[7]) - (1/7)*z[2]**2 - (2/7)*z[2]*z[3] - (2/7)*z[2]*z[4] - (2/7)*z[2]*z[5] - (2/7)*z[2]*z[6] - (2/7)*z[2]*z[7] - (1/9)*z[2]*conj(z[0]) - (1/9)*z[2]*conj(z[1]) - (1/9)*z[2]*conj(z[2]) - (1/9)*z[2]*conj(z[3]) - (1/9)*z[2]*conj(z[4]) - (1/9)*z[2]*conj(z[5]) - (1/9)*z[2]*conj(z[6]) - (1/9)*z[2]*conj(z[7]) - (1/7)*z[3]**2 - (2/7)*z[3]*z[4] - (2/7)*z[3]*z[5] - (2/7)*z[3]*z[6] - (2/7)*z[3]*z[7] - (1/9)*z[3]*conj(z[0]) - (1/9)*z[3]*conj(z[1]) - (1/9)*z[3]*conj(z[2]) - (1/9)*z[3]*conj(z[3]) - (1/9)*z[3]*conj(z[4]) - (1/9)*z[3]*conj(z[5]) - (1/9)*z[3]*conj(z[6]) - (1/9)*z[3]*conj(z[7]) - (1/7)*z[4]**2 - (2/7)*z[4]*z[5] - (2/7)*z[4]*z[6] - (2/7)*z[4]*z[7] - (1/9)*z[4]*conj(z[0]) - (1/9)*z[4]*conj(z[1]) - (1/9)*z[4]*conj(z[2]) - (1/9)*z[4]*conj(z[3]) - (1/9)*z[4]*conj(z[4]) - (1/9)*z[4]*conj(z[5]) - (1/9)*z[4]*conj(z[6]) - (1/9)*z[4]*conj(z[7]) - (1/7)*z[5]**2 - (2/7)*z[5]*z[6] - (2/7)*z[5]*z[7] - (1/9)*z[5]*conj(z[0]) - (1/9)*z[5]*conj(z[1]) - (1/9)*z[5]*conj(z[2]) - (1/9)*z[5]*conj(z[3]) - (1/9)*z[5]*conj(z[4]) - (1/9)*z[5]*conj(z[5]) - (1/9)*z[5]*conj(z[6]) - (1/9)*z[5]*conj(z[7]) - (1/7)*z[6]**2 - (2/7)*z[6]*z[7] - (1/9)*z[6]*conj(z[0]) - (1/9)*z[6]*conj(z[1]) - (1/9)*z[6]*conj(z[2]) - (1/9)*z[6]*conj(z[3]) - (1/9)*z[6]*conj(z[4]) - (1/9)*z[6]*conj(z[5]) - (1/9)*z[6]*conj(z[6]) - (1/9)*z[6]*conj(z[7]) - (1/7)*z[7]**2 - (1/9)*z[7]*conj(z[0]) - (1/9)*z[7]*conj(z[1]) - (1/9)*z[7]*conj(z[2]) - (1/9)*z[7]*conj(z[3]) - (1/9)*z[7]*conj(z[4]) - (1/9)*z[7]*conj(z[5]) - (1/9)*z[7]*conj(z[6]) - (1/9)*z[7]*conj(z[7]) + 1*conj(z[0])**2 - (2/7)*conj(z[0])*conj(z[1]) - (2/7)*conj(z[0])*conj(z[2]) - (2/7)*conj(z[0])*conj(z[3]) - (2/7)*conj(z[0])*conj(z[4]) - (2/7)*conj(z[0])*conj(z[5]) - (2/7)*conj(z[0])*conj(z[6]) - (2/7)*conj(z[0])*conj(z[7]) - (1/7)*conj(z[1])**2 - (2/7)*conj(z[1])*conj(z[2]) - (2/7)*conj(z[1])*conj(z[3]) - (2/7)*conj(z[1])*conj(z[4]) - (2/7)*conj(z[1])*conj(z[5]) - (2/7)*conj(z[1])*conj(z[6]) - (2/7)*conj(z[1])*conj(z[7]) - (1/7)*conj(z[2])**2 - (2/7)*conj(z[2])*conj(z[3]) - (2/7)*conj(z[2])*conj(z[4]) - (2/7)*conj(z[2])*conj(z[5]) - (2/7)*conj(z[2])*conj(z[6]) - (2/7)*conj(z[2])*conj(z[7]) - (1/7)*conj(z[3])**2 - (2/7)*conj(z[3])*conj(z[4]) - (2/7)*conj(z[3])*conj(z[5]) - (2/7)*conj(z[3])*conj(z[6]) - (2/7)*conj(z[3])*conj(z[7]) - (1/7)*conj(z[4])**2 - (2/7)*conj(z[4])*conj(z[5]) - (2/7)*conj(z[4])*conj(z[6]) - (2/7)*conj(z[4])*conj(z[7]) - (1/7)*conj(z[5])**2 - (2/7)*conj(z[5])*conj(z[6]) - (2/7)*conj(z[5])*conj(z[7]) - (1/7)*conj(z[6])**2 - (2/7)*conj(z[6])*conj(z[7]) - (1/7)*conj(z[7])**2
def dbdz(z):
    diff = [2*z[0] - (2/7)*z[1] - (2/7)*z[2] - (2/7)*z[3] - (2/7)*z[4] - (2/7)*z[5] - (2/7)*z[6] - (2/7)*z[7] - 1*conj(z[0]) - (1/9)*conj(z[1]) - (1/9)*conj(z[2]) - (1/9)*conj(z[3]) - (1/9)*conj(z[4]) - (1/9)*conj(z[5]) - (1/9)*conj(z[6]) - (1/9)*conj(z[7])]
    diff += 7*[- (2/7)*z[0] - (2/7)*z[1] - (2/7)*z[2] - (2/7)*z[3] - (2/7)*z[4] - (2/7)*z[5] - (2/7)*z[6] - (2/7)*z[7] - (1/9)*conj(z[0]) - (1/9)*conj(z[1]) - (1/9)*conj(z[2]) - (1/9)*conj(z[3]) - (1/9)*conj(z[4]) - (1/9)*conj(z[5]) - (1/9)*conj(z[6]) - (1/9)*conj(z[7])]
    return diff

def dbdzconj(z):
    diff_conj = [2*conj(z[0]) - (2/7)*conj(z[1]) - (2/7)*conj(z[2]) - (2/7)*conj(z[3]) - (2/7)*conj(z[4]) - (2/7)*conj(z[5]) - (2/7)*conj(z[6]) - (2/7)*conj(z[7])- 1*z[0] - (1/9)*z[1] - (1/9)*z[2] - (1/9)*z[3] - (1/9)*z[4] - (1/9)*z[5] - (1/9)*z[6] - (1/9)*z[7]]
    diff_conj += 7*[- (2/7)*conj(z[0]) - (2/7)*conj(z[1]) - (2/7)*conj(z[2]) - (2/7)*conj(z[3]) - (2/7)*conj(z[4]) - (2/7)*conj(z[5]) - (2/7)*conj(z[6]) - (2/7)*conj(z[7]) - (1/9)*z[0] - (1/9)*z[1] - (1/9)*z[2] - (1/9)*z[3] - (1/9)*z[4] - (1/9)*z[5] - (1/9)*z[6] - (1/9)*z[7]]
    return diff_conj

# Defining barrier and Hamiltonian for n
b = barrier
ham = H_hat_n(n)

# Conditions for barrier
init = lambda z: Implies(norm_sqrd(z[0]) >= 0.9, And(b(z).r <= 0, b(z).i == 0))
# unsafe = lambda z: Implies(Or([norm_sqrd(zi) >= 0.9 for zi in z[1:]]), And(b(z).r > 0, b(z).i == 0))
unsafe = lambda z: Implies(norm_sqrd(z[1]) >= 0.9, And(b(z).r > 0, b(z).i == 0))
dbdt = lambda z: np.dot(dbdz(z), evo(ham, z)) + np.dot(dbdzconj(z), [conj(e) for e in evo(ham, z)])
evo_cond = lambda z: simplify(dbdt(z).r <= 0)
evo_isreal_cond = lambda z: simplify(dbdt(z).i == 0)
conds = [
    evo_isreal_cond,
    evo_cond,
    init,
    unsafe,
    ]
strs = [
    "Evolution is real",
    "Evolution",
    "Init",
    "Unsafe",
    ]

if __name__ == "__main__":
    s = z3.Solver()
    z = ComplexVector('z', N)
    s.add(sum_norm_sqrs(z) == 1)
    s.push()
    for cond, cond_str in zip(conds, strs):
        s.add(Not(cond(z)))
        print("Checking " + cond_str + " condition...")
        sat = s.check()

        if sat == z3.sat:
            print("Found CEx")
            mz = [0] * N
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
            exit()
        else: print(sat)

        print()
        s.pop()