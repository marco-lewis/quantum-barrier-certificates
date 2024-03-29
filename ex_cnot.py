from src.find_barrier import *
from src.utils import *

H = [[0,0,0,0],[0,0,0,0],[0,0,1,-1],[0,0,-1,1]]
k = 2

print("5.3: 1a. Start near 00")
def init(x): return [x[0]**2 + x[1]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
def unsafe(x): return [x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.11,1], [1,1])]

expr = z0*z0.conjugate()
b_cnot = scipy_find_k_barrier(k, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              verbose=1,
                              objective_expressions=[expr]
                              )

print("5.3: 1b. Start near 01")
def init(x): return [x[2]**2 + x[3]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0.9, 1], [1, 1])]
def unsafe(x): return [x[0]**2 + x[1]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.11,1], [1,1])]

expr = z1*z1.conjugate()
b_cnot = scipy_find_k_barrier(k, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              verbose=1,
                              objective_expressions=[expr]
                              )

print("5.3: 2a. Start near 10")
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

expr = z3*z3.conjugate()
b_cnot = scipy_find_k_barrier(k, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              verbose=1,
                              objective_expressions=[expr]
                              )

print("5.3: 2b. Start near 11")
def init(x): return [x[0]**2 + x[1]**2,
                     x[2]**2 + x[3]**2,
                     x[4]**2 + x[5]**2,
                     x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
init_constraints = [NonlinearConstraint(init, [0, 0, 0, 0.81, 1], [0.01, 0.09, 0.09, 1, 1])]
def unsafe(x): return [x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
                     x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2]
unsafe_constraints= [NonlinearConstraint(unsafe, [.5,0,1], [1,.5,1])]

b_cnot = scipy_find_k_barrier(k, H,
                              init=init_constraints,
                              unsafe=unsafe_constraints,
                              verbose=1,
                              objective_expressions=[expr]
                              )