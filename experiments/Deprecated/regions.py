import numpy as np

# 1-qubit regions
def near_0():
    r = np.random.uniform(np.sqrt(0.9), 1)
    ang0 = np.random.uniform(0, 2*np.pi)
    ang1 = np.random.uniform(0, 2*np.pi)
    z0 = r * np.exp(1j*ang0)
    z1 = np.sqrt(1 - r**2) * np.exp(1j*ang1)
    return [z0, z1]

def near_1():
    z = near_0()
    return [z[1], z[0]]

def near_0_sup():
    r = np.random.uniform(np.sqrt(0.9), 1)
    ang0 = np.random.uniform(0, 2*np.pi)
    ang1 = np.random.uniform(0, 2*np.pi)
    z0 = r * np.exp(1j*ang0)
    z1 = np.sqrt(1 - r**2) * np.exp(1j*ang1)
    return [z0, z1]

def near_1_suped():
    z = near_0_sup()
    return [z[1], z[0]]

# Phased regions
a = .2
def near_unphased():
    r = np.random.uniform(0, 1)
    ang1 = np.random.uniform(-a, a)
    delta = np.random.uniform(0, 2*np.pi)
    z0 = r * np.exp(1j*delta)
    z1 = np.sqrt(1 - r**2) * np.exp(1j*(ang1 + delta))
    return [z0, z1]

def very_phased():
    r = np.random.uniform(0, 1)
    ang1 = np.random.uniform(np.pi - a, np.pi + a)
    delta = np.random.uniform(0, 2*np.pi)
    z0 = r * np.exp(1j*delta)
    z1 = np.sqrt(1 - r**2) * np.exp(1j*(ang1 + delta))
    return [z0, z1]

# 2-qubit regions
def ctrl_near_0():
    z = near_0()
    return z[:1] + [0] + z[1:] + [0]

def tgt_near_0():
    z = near_0()
    return z[:1] + z[1:] + [0]*2

def near_00():
    ctrl, tgt = near_0(), near_0()
    z = []
    for c in ctrl:
        for t in tgt:
            z.append(c*t)
    return z

def near_01():
    ctrl, tgt = near_0(), near_0()
    z = [ctrl[0] * tgt[1], ctrl[0]*tgt[0], ctrl[1] * tgt[1], ctrl[1] * tgt[0]]
    return z

def near_10():
    ctrl, tgt = near_0(), near_0()
    z = [ctrl[1] * tgt[0], ctrl[1]*tgt[1], ctrl[0] * tgt[0], ctrl[0] * tgt[1]]
    return z

def away_from_11_10():
    if np.random.rand() > 0.5: return near_00()
    else: return near_01()

def no_reach_00_01():
    pass
    
def near_11():
    r = np.random.uniform(np.sqrt(0.81),1)
    not_r = np.sqrt(1 - r**2)
    nr = np.random.uniform(0,1, (3,))
    nr *= not_r/(nr**2).sum()
    rs = nr.tolist()
    rs.insert(3,r)
    angs = np.random.uniform(0, 2*np.pi, (4,)).tolist()
    out = []
    for r, a in zip(rs, angs):
        out.append(r * np.exp(1j*a))
    return out
