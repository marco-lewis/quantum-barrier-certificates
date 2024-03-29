from utils_plot import *

def cartesian_evo(t, x0):
    x, y, z = x0
    mid = (x + z)/ 2
    dist = (z - x) /2
    cost = np.cos(2*t)
    sint = np.sin(2*t)
    
    # Define functions
    x = mid - dist * cost
    y = - np.sqrt(2) * dist * sint
    z = mid + dist * cost
    return (x,y,z)


# # Hadamard Hamiltonian
def cart_bloch_near_0():
    phi, delta = randang(), randang()
    r = np.random.uniform(np.sqrt(0.9), 1)
    # State form
    z0 = np.array([r * np.exp(1j * delta), np.sqrt(1 - r**2) * np.exp(1j * phi)])
    x0 = state_to_bloch([z0[0]], [z0[1]])
    # Cartesian form
    X0 = x0[0] + x0[1] + x0[2]
    return X0, z0

def cart_bloch_near_1():
    phi, delta = randang(), randang()
    r = np.random.uniform(0, np.sqrt(0.1))
    # State form
    z0 = np.array([r * np.exp(1j * delta), np.sqrt(1 - r**2) * np.exp(1j * phi)])
    x0 = state_to_bloch([z0[0]], [z0[1]])
    # Cartesian form
    X0 = x0[0] + x0[1] + x0[2]
    return X0, z0

def start_near_0():
    # State form
    z0 = np.array([np.sqrt(0.9), -np.sqrt(0.1)])
    x0 = state_to_bloch([z0[0]], [z0[1]])
    # Cartesian form
    X0 = x0[0] + x0[1] + x0[2]
    return X0, z0

def start_near_1():
    # State form
    z0 = np.array([np.sqrt(0.1), -np.sqrt(0.9)])
    x0 = state_to_bloch([z0[0]], [z0[1]])
    # Cartesian form
    X0 = x0[0] + x0[1] + x0[2]
    return X0, z0

start_points_0 = [([-0.4205748172100388, -0.19077787500541843, 0.8869727309997525], np.array([-0.57843339+0.7803212j ,  0.20781569-0.11543948j]))]
start_points_1 = [([-0.019089171625820144, -0.46005262894314497, -0.8876864210570876], np.array([-0.03054791+0.23499705j,  0.96777386+0.08518778j]))]

step = 0.1
t0 = 0
t = 2*np.pi
dXdt = lambda x, t: np.array([-x[1], (x[0] - x[2]), x[1]])

trange = np.arange(t0, t, step=step)
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111,projection='3d',box_aspect=(1,1,1))

barrier=True
plot_H_bloch_regions(ax, barrier)

for X0, Z0 in start_points_0:
    dataSet, ts = rk4(X0, t0, t, step, dXdt)
    dataSet = dataSet.T
    xps, yps, zps = dataSet[0], dataSet[1], dataSet[2]

    plt.quiver(xps[:-1], yps[:-1], zps[:-1], xps[1:]-xps[:-1], yps[1:]-yps[:-1], zps[1:]-zps[:-1], color='black')
    ax.scatter(dataSet.T[0][0],dataSet.T[0][1],dataSet.T[0][2], label=str(dataSet.T[0]), s=100, c='black')
    a = cartesian_evo(np.pi/2,X0)
    b = cartesian_evo(np.pi/2 + 0.01,X0)
    
    z_vec, o_vec = get_states(trange, Z0, had_ham)
    xs, ys, zs, _ = state_to_bloch(z_vec, o_vec)


fig.tight_layout()

if not barrier: save("had_no_barrier.png")
save('2a_had_bloch.png')
ax.set_yticks([])
ax.set_ylabel('')
ax.azim = -90
ax.elev = 1
save('2b_had_side_bloch.png')

ax.set_zticks([])
ax.set_zlabel('')
ax.azim = 270
ax.elev = 90
# save('had_above_bloch.png')
plt.close()
print("Hadamard figures stored")

# # S-gate Hamiltonian
step = 0.1
t0 = 0
t = 2*np.pi
dXdt = lambda x, t: np.array([-x[1], x[0], 0])

trange = np.arange(t0, t, step=step)
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111,projection='3d',box_aspect=(1,1,1))

plot_S_bloch_regions(ax, True)

for X0, Z0 in start_points_0:
    dataSet, ts = rk4(X0, t0, t, step, dXdt)
    dataSet = dataSet.T
    xps, yps, zps = dataSet[0], dataSet[1], dataSet[2]

    # ax.plot(xps,yps,zps, label=str(X0), color='black', ls='--')
    plt.quiver(xps[:-1], yps[:-1], zps[:-1], xps[1:]-xps[:-1], yps[1:]-yps[:-1], zps[1:]-zps[:-1], color='black')
    ax.scatter(dataSet.T[0][0],dataSet.T[0][1],dataSet.T[0][2], label=str(dataSet.T[0]), s=100, c='black')
    a = cartesian_evo(np.pi/2,X0)
    b = cartesian_evo(np.pi/2 + 0.01,X0)
    
    z_vec, o_vec = get_states(trange, Z0, had_ham)
    xs, ys, zs, _ = state_to_bloch(z_vec, o_vec)

fig.tight_layout()
# save('S0_normal_bloch.png')

ax.set_yticks([])
ax.set_ylabel('')
ax.azim = -90
ax.elev = 2
save('3a_S0_side_bloch.png')

ax.set_zticks([])
ax.set_zlabel('')
ax.azim = 270
ax.elev = 90
# save('S0_above_bloch.png')
plt.close()
print("Phase on |0> figures stored")

step = 0.1
t0 = 0
t = 2*np.pi
dXdt = lambda x, t: np.array([-x[1], x[0], 0])

trange = np.arange(t0, t, step=step)
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111,projection='3d',box_aspect=(1,1,1))

plot_S_bloch_regions(ax, True, init=1)

for X0, Z0 in start_points_1:
    dataSet, ts = rk4(X0, t0, t, step, dXdt)
    dataSet = dataSet.T
    xps, yps, zps = dataSet[0], dataSet[1], dataSet[2]

    plt.quiver(xps[:-1], yps[:-1], zps[:-1], xps[1:]-xps[:-1], yps[1:]-yps[:-1], zps[1:]-zps[:-1], color='black')
    ax.scatter(dataSet.T[0][0],dataSet.T[0][1],dataSet.T[0][2], label=str(dataSet.T[0]), s=100, c='black')
    a = cartesian_evo(np.pi/2,X0)
    b = cartesian_evo(np.pi/2 + 0.01,X0)
    
    z_vec, o_vec = get_states(trange, Z0, had_ham)
    xs, ys, zs, _ = state_to_bloch(z_vec, o_vec)


fig.tight_layout()

# save('S1_normal_bloch.png')
ax.set_yticks([])
ax.set_ylabel('')
ax.azim = -90
ax.elev = -2
save('3b_S1_side_bloch.png')

ax.set_zticks([])
ax.set_zlabel('')
ax.azim = 270
ax.elev = 90
# save('S1_above_bloch.png')

plt.close()
print("Phase on |1> figures stored")