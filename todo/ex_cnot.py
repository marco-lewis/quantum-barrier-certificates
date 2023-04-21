from src.find_barrier import *
from src.utils import *

print("4. Slightly Disturbed")
d = 0.01
H = [[0,0,0,0],
     [0,0,0,0],
     [0,0,1,-1],
     [0,0,-1,1]]
i, j = 1,2
H[i][j] = d
H[j][i] = d


def init(x): return [x[0]**2 + x[1]**2,
                     x[2]**2 + x[3]**2,
                     x[4]**2 + x[5]**2,
                     x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0, 0, 0.81, 0, 1], [0.09, 0.01, 1, 0.09, 1])]
def unsafe(x): return [x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
                     x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.5,0,1], [1,.5,1])]
try:
    b_ctop = scipy_find_k_barrier(2, H,
                                  init=init_constraints,
                                  unsafe=unsafe_constraints,
                                  verbose=1)
    print(b_ctop)
except: print("Fail")