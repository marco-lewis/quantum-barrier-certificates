from find_barrier import *
from utils import *

n = 2
H_had = [[1,1],[1,-1]]
kron_had = lambda H: np.kron(H, H_had)
H = H_had
for i in range(n-1): H = kron_had(H)

sos = lambda xs: sum([x**2 for x in xs])

def init(x): return [x[0]**2 + x[1]**2, sos(x)]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]

def unsafe(x): return [x[2*i]**2 + x[2*i + 1]**2 for i in range(1,2**n)] + [sos(x)]
unsafe_constraints= [NonlinearConstraint(unsafe, ([.5] * (2**n - 1)) + [1], [1] * 2**n)]

coeff_num = len(generate_term_powers(2, 2))
bs = []
try:
    b_ctop = scipy_find_k_barrier(2, H,
                                      init=init_constraints,
                                      unsafe=unsafe_constraints,
                                      # linprog_obj=[],
                                      verbose=1)
    bs.append(b_ctop)
except:
    print("Fail")
print(bs)