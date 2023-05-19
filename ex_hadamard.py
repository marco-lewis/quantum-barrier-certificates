from src.find_barrier import *
from src.utils import *

# See Section 4.1 of https://doi.org/10.1088/1361-6455/aa987c
# Hadamard as a Hamiltonian
# H = [[1/2*np.sqrt(2), 1/2*np.sqrt(2)],[1/2*np.sqrt(2), -1/2*np.sqrt(2)]]
H = [[1,1],[1,-1]]
k = 2

print("5.1: Hadamard Barrier")
init_constraints = close_to_zero
unsafe_constraints = qubit_simple_constraints((0,0.1),(0.9,1))

barrier = scipy_find_k_barrier(k, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1,
                               obj_func_idxs=[12]
                               )

barrier = 3*barrier
print("Scaled: ", barrier)
z0 = sym.Symbol("z0", complex=True)
z1 = sym.Symbol("z1", complex=True)
expr = 1 - z0*z0.conjugate()
print("Substituting z1*conj(z1) = 1 - z0*conj(z0): ", barrier.subs(z1*z1.conjugate(), expr))