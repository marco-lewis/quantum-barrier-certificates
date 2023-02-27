from find_barrier import *
from utils import *

# See https://doi.org/10.1088/1361-6455/aa987c
# Hadamard as a Hamiltonian
# H = [[1/2*np.sqrt(2), 1/2*np.sqrt(2)],[1/2*np.sqrt(2), -1/2*np.sqrt(2)]]
H = [[1,1],[1,-1]]
k = 2

init_constraints = close_to_zero
unsafe_constraints = qubit_simple_constraints((0,0.1),(0.9,1))

barrier = scipy_find_k_barrier(k, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1)

print_barrier(barrier)
print("Scaled: ", 3*barrier)