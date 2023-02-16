from general_had import *

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
sum_term = lambda z: sum([z[0] * conj(zi) for zi in z[1:]]) + sum([conj(z[0]) * zi for zi in z[1:]])
barrier = lambda z: 1.7 - 2 * z[0] * conj(z[0]) - sum_term(z)

def dbdz(z):
    diff = [-2 * conj(z[0]) - np.sum([conj(zi) for zi in z[1:]])]
    for j in range(1, len(z)): diff.append(0 - conj(z[0]))
    return diff

def dbdzconj(z):
    diff_conj = [-2 * z[0] - np.sum(z[1:])]
    for j in range(1, len(z)): diff_conj.append(0 - z[0])
    return diff_conj

# Attempt 3

# Attempt n