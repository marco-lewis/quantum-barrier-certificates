from find_barrier import *
from utils import *

print("1. Start near 0")
H = [[0,0,0,0],[0,0,0,0],[0,0,1,-1],[0,0,-1,1]]
def init(x): return [x[0]**2 + x[1]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
def unsafe(x): return [x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.11,1], [1,1])]

b_cnot = scipy_find_k_barrier(2, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              linprog_obj=[40],
                              verbose=1)
print_barrier(b_cnot)

print("2. Start near 00 or 01")
def init(x): return [x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
def unsafe(x): return [x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.11,1], [1,1])]

b_cnot = scipy_find_k_barrier(2, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              linprog_obj=[40, 31],
                              verbose=1)
print_barrier(b_cnot)

print("3a. Start near 01")
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

b_cnot = scipy_find_k_barrier(2, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              linprog_obj=[16],
                              verbose=1)
print_barrier(b_cnot)

print("3b. Start near 01 (slightly different init constraints)")
def init(x): return [x[4]**2 + x[5]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
def unsafe(x): return [x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.11,1], [1,1])]

b_cnot = scipy_find_k_barrier(2, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              linprog_obj=[16],
                              verbose=1)
print_barrier(b_cnot)

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
                                  linprog_obj=[16,23],
                                #   linprog_obj=[71,130],
                                  verbose=1)
    print(b_ctop)
except: print("Fail")