from todo.general_had import *

# Attempt 1
term1 = lambda z: np.prod([conj(z[i]) if i > 0 else z[i] for i in range(len(z))])
term2 = lambda z: np.prod([z[i] if i > 0 else conj(z[i]) for i in range(len(z))])
barrier = lambda z: 6/5 - 2 * z[0] * conj(z[0]) - term1(z) - term2(z)

def dbdz(z):
    diff = [-2 * conj(z[0]) - np.prod([conj(zi) for zi in z][1:])]
    for j in range(1, len(z)): diff.append(0 - np.prod([1 if i == j else z[i] if i > 0 else conj(z[i]) for i in range(len(z))]))
    return diff

def dbdzconj(z):
    diff_conj = [-2 * z[0] - np.prod(z[1:])]
    for j in range(1, len(z)): diff_conj.append(0 - np.prod([1 if i == j else conj(z[i]) if i > 0 else z[i] for i in range(len(z))]))
    return diff_conj

# Attempt 2
z_coeff = 2
sum_term = lambda z: sum([z[0] * conj(zi) for zi in z[1:]]) + sum([conj(z[0]) * zi for zi in z[1:]])
barrier = lambda z: 1.2 - (z_coeff * z[0] * conj(z[0]) + sum_term(z))

def dbdz(z):
    diff = [-z_coeff * conj(z[0]) - np.sum([conj(zi) for zi in z[1:]])]
    for j in range(1, len(z)): diff.append(0 - conj(z[0]))
    return diff

def dbdzconj(z):
    diff_conj = [-z_coeff * z[0] - np.sum(z[1:])]
    for j in range(1, len(z)): diff_conj.append(0 - z[0])
    return diff_conj

# Attempt 3
sum_term = lambda z0, z1: z0 * conj(z1) + conj(z0) * z1
barrier = lambda z: 1.2 - sum([sum([sum_term(zi, zj) for zj in z]) for zi in z])
# Differential is prop. to sum(-Im(z_0 * z_x) for zx in z)

# Attempt 4
# Generated from Scipy
# Requires n = 3
c = -1
barrier = lambda z: c + 1*z[0]**2 - (2/7)*z[0]*z[1] - (2/7)*z[0]*z[2] - (2/7)*z[0]*z[3] - (2/7)*z[0]*z[4] - (2/7)*z[0]*z[5] - (2/7)*z[0]*z[6] - (2/7)*z[0]*z[7] - 1*z[0]*conj(z[0]) - (1/9)*z[0]*conj(z[1]) - (1/9)*z[0]*conj(z[2]) - (1/9)*z[0]*conj(z[3]) - (1/9)*z[0]*conj(z[4]) - (1/9)*z[0]*conj(z[5]) - (1/9)*z[0]*conj(z[6]) - (1/9)*z[0]*conj(z[7]) - (1/7)*z[1]**2 - (2/7)*z[1]*z[2] - (2/7)*z[1]*z[3] - (2/7)*z[1]*z[4] - (2/7)*z[1]*z[5] - (2/7)*z[1]*z[6] - (2/7)*z[1]*z[7] - (1/9)*z[1]*conj(z[0]) - (1/9)*z[1]*conj(z[1]) - (1/9)*z[1]*conj(z[2]) - (1/9)*z[1]*conj(z[3]) - (1/9)*z[1]*conj(z[4]) - (1/9)*z[1]*conj(z[5]) - (1/9)*z[1]*conj(z[6]) - (1/9)*z[1]*conj(z[7]) - (1/7)*z[2]**2 - (2/7)*z[2]*z[3] - (2/7)*z[2]*z[4] - (2/7)*z[2]*z[5] - (2/7)*z[2]*z[6] - (2/7)*z[2]*z[7] - (1/9)*z[2]*conj(z[0]) - (1/9)*z[2]*conj(z[1]) - (1/9)*z[2]*conj(z[2]) - (1/9)*z[2]*conj(z[3]) - (1/9)*z[2]*conj(z[4]) - (1/9)*z[2]*conj(z[5]) - (1/9)*z[2]*conj(z[6]) - (1/9)*z[2]*conj(z[7]) - (1/7)*z[3]**2 - (2/7)*z[3]*z[4] - (2/7)*z[3]*z[5] - (2/7)*z[3]*z[6] - (2/7)*z[3]*z[7] - (1/9)*z[3]*conj(z[0]) - (1/9)*z[3]*conj(z[1]) - (1/9)*z[3]*conj(z[2]) - (1/9)*z[3]*conj(z[3]) - (1/9)*z[3]*conj(z[4]) - (1/9)*z[3]*conj(z[5]) - (1/9)*z[3]*conj(z[6]) - (1/9)*z[3]*conj(z[7]) - (1/7)*z[4]**2 - (2/7)*z[4]*z[5] - (2/7)*z[4]*z[6] - (2/7)*z[4]*z[7] - (1/9)*z[4]*conj(z[0]) - (1/9)*z[4]*conj(z[1]) - (1/9)*z[4]*conj(z[2]) - (1/9)*z[4]*conj(z[3]) - (1/9)*z[4]*conj(z[4]) - (1/9)*z[4]*conj(z[5]) - (1/9)*z[4]*conj(z[6]) - (1/9)*z[4]*conj(z[7]) - (1/7)*z[5]**2 - (2/7)*z[5]*z[6] - (2/7)*z[5]*z[7] - (1/9)*z[5]*conj(z[0]) - (1/9)*z[5]*conj(z[1]) - (1/9)*z[5]*conj(z[2]) - (1/9)*z[5]*conj(z[3]) - (1/9)*z[5]*conj(z[4]) - (1/9)*z[5]*conj(z[5]) - (1/9)*z[5]*conj(z[6]) - (1/9)*z[5]*conj(z[7]) - (1/7)*z[6]**2 - (2/7)*z[6]*z[7] - (1/9)*z[6]*conj(z[0]) - (1/9)*z[6]*conj(z[1]) - (1/9)*z[6]*conj(z[2]) - (1/9)*z[6]*conj(z[3]) - (1/9)*z[6]*conj(z[4]) - (1/9)*z[6]*conj(z[5]) - (1/9)*z[6]*conj(z[6]) - (1/9)*z[6]*conj(z[7]) - (1/7)*z[7]**2 - (1/9)*z[7]*conj(z[0]) - (1/9)*z[7]*conj(z[1]) - (1/9)*z[7]*conj(z[2]) - (1/9)*z[7]*conj(z[3]) - (1/9)*z[7]*conj(z[4]) - (1/9)*z[7]*conj(z[5]) - (1/9)*z[7]*conj(z[6]) - (1/9)*z[7]*conj(z[7]) + 1*conj(z[0])**2 - (2/7)*conj(z[0])*conj(z[1]) - (2/7)*conj(z[0])*conj(z[2]) - (2/7)*conj(z[0])*conj(z[3]) - (2/7)*conj(z[0])*conj(z[4]) - (2/7)*conj(z[0])*conj(z[5]) - (2/7)*conj(z[0])*conj(z[6]) - (2/7)*conj(z[0])*conj(z[7]) - (1/7)*conj(z[1])**2 - (2/7)*conj(z[1])*conj(z[2]) - (2/7)*conj(z[1])*conj(z[3]) - (2/7)*conj(z[1])*conj(z[4]) - (2/7)*conj(z[1])*conj(z[5]) - (2/7)*conj(z[1])*conj(z[6]) - (2/7)*conj(z[1])*conj(z[7]) - (1/7)*conj(z[2])**2 - (2/7)*conj(z[2])*conj(z[3]) - (2/7)*conj(z[2])*conj(z[4]) - (2/7)*conj(z[2])*conj(z[5]) - (2/7)*conj(z[2])*conj(z[6]) - (2/7)*conj(z[2])*conj(z[7]) - (1/7)*conj(z[3])**2 - (2/7)*conj(z[3])*conj(z[4]) - (2/7)*conj(z[3])*conj(z[5]) - (2/7)*conj(z[3])*conj(z[6]) - (2/7)*conj(z[3])*conj(z[7]) - (1/7)*conj(z[4])**2 - (2/7)*conj(z[4])*conj(z[5]) - (2/7)*conj(z[4])*conj(z[6]) - (2/7)*conj(z[4])*conj(z[7]) - (1/7)*conj(z[5])**2 - (2/7)*conj(z[5])*conj(z[6]) - (2/7)*conj(z[5])*conj(z[7]) - (1/7)*conj(z[6])**2 - (2/7)*conj(z[6])*conj(z[7]) - (1/7)*conj(z[7])**2
def dbdz(z):
    diff = [2*z[0] - (2/7)*z[1] - (2/7)*z[2] - (2/7)*z[3] - (2/7)*z[4] - (2/7)*z[5] - (2/7)*z[6] - (2/7)*z[7] - 1*conj(z[0]) - (1/9)*conj(z[1]) - (1/9)*conj(z[2]) - (1/9)*conj(z[3]) - (1/9)*conj(z[4]) - (1/9)*conj(z[5]) - (1/9)*conj(z[6]) - (1/9)*conj(z[7])]
    diff += 7*[- (2/7)*z[0] - (2/7)*z[1] - (2/7)*z[2] - (2/7)*z[3] - (2/7)*z[4] - (2/7)*z[5] - (2/7)*z[6] - (2/7)*z[7] - (1/9)*conj(z[0]) - (1/9)*conj(z[1]) - (1/9)*conj(z[2]) - (1/9)*conj(z[3]) - (1/9)*conj(z[4]) - (1/9)*conj(z[5]) - (1/9)*conj(z[6]) - (1/9)*conj(z[7])]
    return diff

def dbdzconj(z):
    diff_conj = [2*conj(z[0]) - (2/7)*conj(z[1]) - (2/7)*conj(z[2]) - (2/7)*conj(z[3]) - (2/7)*conj(z[4]) - (2/7)*conj(z[5]) - (2/7)*conj(z[6]) - (2/7)*conj(z[7])- 1*z[0] - (1/9)*z[1] - (1/9)*z[2] - (1/9)*z[3] - (1/9)*z[4] - (1/9)*z[5] - (1/9)*z[6] - (1/9)*z[7]]
    diff_conj += 7*[- (2/7)*conj(z[0]) - (2/7)*conj(z[1]) - (2/7)*conj(z[2]) - (2/7)*conj(z[3]) - (2/7)*conj(z[4]) - (2/7)*conj(z[5]) - (2/7)*conj(z[6]) - (2/7)*conj(z[7]) - (1/9)*z[0] - (1/9)*z[1] - (1/9)*z[2] - (1/9)*z[3] - (1/9)*z[4] - (1/9)*z[5] - (1/9)*z[6] - (1/9)*z[7]]
    return diff_conj

# Attempt n
