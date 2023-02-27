from math import comb

from find_barrier import *
from utils import *

print("1. Hadamard-n start near 0")
n = 3
k = 2
N = 2**n
H_had = [[1,1],[1,-1]]

print("Making Operation")
kron_had = lambda H: np.kron(H, H_had)
H = H_had
for i in range(n-1): H = kron_had(H)

sos = lambda xs: sum([x**2 for x in xs])

print("Making Constraints")
def init(x): return [x[0]**2 + x[1]**2, sos(x)]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
# All basis states unsafe
# def unsafe(x): return [x[2*i]**2 + x[2*i + 1]**2 for i in range(1,2**n)] + [sos(x)]
# unsafe_constraints= [NonlinearConstraint(unsafe, ([.9] * (2**n - 1)) + [1], [1] * 2**n)]

# |q> state unsafe
q = 1
def unsafe(x): return [x[2*q]**2 + x[2*q + 1]**2] + [sos(x)]
unsafe_constraints= [NonlinearConstraint(unsafe, [.9] + [1], [1] * 2)]

coeff_num = 1
for i in range(1, k + 1): coeff_num += comb(k, i) * comb(2*N, i)

if 0:
    tp = generate_term_powers(k, 2**n)
    for i in range(coeff_num): print(i, tp[i])

bs = []
print("Finding barrier")
try:
    b_ctop = scipy_find_k_barrier(2, H,
                                      init=init_constraints,
                                      unsafe=unsafe_constraints,
                                      verbose=1)
    bs.append(b_ctop)
except:
    print("Fail")
print(bs)