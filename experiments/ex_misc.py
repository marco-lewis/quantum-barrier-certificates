from find_barrier import *
from utils import *

print("1. Hamiltonian for X")
H = [[.5,-.5],[-.5,.5]]
# Need suitable init and unsafe regions
def init(x): return [x[0]**2,
                     x[1]**2,
                     x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2]
init_constraints = [NonlinearConstraint(init, [.9, 0, 1], [1, 0, 1])]

def unsafe(x): return [x[1],
                       x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2]
unsafe_constraints = [NonlinearConstraint(unsafe, [.6, 1], [1, 1])]

try:
    barrier = scipy_find_k_barrier(2, H,
                                   init=init_constraints,
                                   unsafe=unsafe_constraints,
                                   linprog_obj=[0,1],
                                   # linprog_obj=list(range(14)),
                                   verbose=1)
    print_barrier(barrier)
except:
    print("Barrier was in unsafe region")

print("2. Hamiltonian for Y")
H = [[.5, .5j],[-.5j,.5]]
init_constraints = close_to_zero
unsafe_constraints = qubit_simple_constraints((0,1),(0,1))
print("Incomplete\n")

print("3. Hamiltonian for Z")
print("a. Init close to 0")
H = [[2,0],[0,1]]
init_constraints = close_to_zero
unsafe_constraints = qubit_simple_constraints((0,0.89),(0.11,1))

barrier = scipy_find_k_barrier(2, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               linprog_obj=[0,12],
                               verbose=1)
print_barrier(barrier)

print("b. Init close to 1")
init_constraints = close_to_one
unsafe_constraints = qubit_simple_constraints((0.11,1),(0,0.89))

barrier = scipy_find_k_barrier(2, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               linprog_obj=[0,7],
                               verbose=1)
barrier