from src.find_barrier import *
from src.utils import *

d=1
print("3. Slight Disturbance")
H = [[d,0.1], [0.1,-d]]
init_constraints = close_to_one
unsafe_constraints = qubit_simple_constraints((np.sqrt(0.11),1),(0,np.sqrt(0.89)))

barrier = scipy_find_k_barrier(2, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1)
print_barrier(barrier)

print("4. Bounded Evolution Experiment")
H = [[d,0], [0,-d]]
print("a. Barrier Above")
init_constraints = qubit_simple_constraints((np.sqrt(0.4),np.sqrt(0.6)),(np.sqrt(0.4),np.sqrt(0.6)))
unsafe_constraints = qubit_simple_constraints((np.sqrt(0.7),1),(0,np.sqrt(0.3)))
barrier = scipy_find_k_barrier(2, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1)
print_barrier(barrier)

print("b. Barrier Below")
init_constraints = qubit_simple_constraints((np.sqrt(0.4),np.sqrt(0.6)),(np.sqrt(0.4),np.sqrt(0.6)))
unsafe_constraints = qubit_simple_constraints((0,np.sqrt(0.3)),(np.sqrt(0.7),1))
barrier = scipy_find_k_barrier(2, H,
                               init=init_constraints,
                               unsafe=unsafe_constraints,
                               verbose=1)
print_barrier(barrier)