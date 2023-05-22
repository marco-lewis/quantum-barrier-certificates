from src.find_barrier import *
from src.utils import *

# Pauli-Z gate as a Hamiltonian
H = [[1,0], [0,-1]]
k = 2

print("5.2: 1. Close to Zero experiment")
init_constraints = close_to_zero
unsafe_constraints = qubit_simple_constraints((0,np.sqrt(0.89)),(np.sqrt(0.11),1))

barrier = scipy_find_k_barrier(k, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1,
                               obj_func_idxs=[12]
                               )

print("5.2: 2. Close to One experiment")
init_constraints = close_to_one
unsafe_constraints = qubit_simple_constraints((np.sqrt(0.11),1),(0,np.sqrt(0.89)))

barrier = scipy_find_k_barrier(k, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1,
                               obj_func_idxs=[7]
                               )