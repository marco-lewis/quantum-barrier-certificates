import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.linalg import expm

import os

import cmasher as cmr
from IPython.display import HTML
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

trange = np.arange(0, 2*np.pi, 0.01)

isqrt2 = 1/np.sqrt(2)
had = np.array([[isqrt2, isqrt2], [isqrt2, -isqrt2]])
X = np.array([[0,1],[1,0]])
Y = np.array([[0,1j],[-1j,0]])
Z = np.array([[1,0],[0,-1]])

def ham(op):
    return lambda t: expm(-1j *op * t)

def had_ham(t):
    return expm(-1j * had * t)

prop = mpl.font_manager.FontProperties()

path = os.path.dirname(os.path.realpath(__file__))
save = lambda file: plt.savefig(path + "/Images/" + file, bbox_inches='tight', pad_inches=0)

def plt_transformation(trange, x_0, ham, ax, ax2):
    xs = [np.dot(ham(t), x_0) for t in trange]
    data = [x[0] for x in xs]
    # extract real part
    x0 = [ele.real for ele in data]
    # extract imaginary part
    y0 = [ele.imag for ele in data]

    data = [x[1] for x in xs]
    x1 = [ele.real for ele in data]
    y1 = [ele.imag for ele in data]
    
    ax.plot(x0,y0, label=str(x_0))
    ax2.plot(x1,y1, linestyle='--', label=str(x_0))
    ax.plot(x_0[0].real, x_0[0].imag, 'x', color='black', markersize=20)
    ax2.plot(x_0[1].real, x_0[1].imag, 'x', color='black', markersize=20)

def randang(): return np.random.uniform(0, 2*np.pi)

def plot_H_bloch_regions(ax, barrier=False):
    state_size = 300
    # Plot |0> and |1>
    ax.scatter(0,0,1.05, label='|0>', marker='^', s=state_size, c='black', zorder=100)
    ax.scatter(0,0,-1.05, label='|1>', marker='v', s=state_size, c='black', zorder=100)
    
    # Plot sphere
    phi = np.linspace(0, np.pi, 20)
    theta = np.linspace(0, 2 * np.pi, 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='w', alpha=0.3, linewidth=0)
    
    # Plot initial
    phi = np.linspace(0, 2*np.pi, 20)
    theta = np.linspace(0, 2*np.arccos(np.sqrt(0.9)), 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x,y,z, color='g', alpha=0.7)
    
    # Plot unsafe
    theta = np.linspace(np.pi, 2*np.arccos(np.sqrt(0.1)), 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x,y,z, color='r', alpha=0.7)
    
    # Plot barrier
    if barrier:
        X,Y = np.meshgrid(np.arange(-1.1,1.1,step=0.1), np.arange(-1.1,1.1,step=0.1))
        f = 0
        # Edge cases (closest to unsafe/initial resp.)
        # f = .4**2
        f = -(.4**2)
        Z = np.array([[f] * len(X[0])] * len(X))
        xp = 1/np.sqrt(2) * (X -Z)
        yp = Y
        zp = 1/np.sqrt(2) * -(X + Z)
        ax.plot_surface(xp,yp,zp, alpha=0.2)

    ax.set_xlabel(u'x', fontproperties=prop, rotation=0, labelpad=10)
    ax.set_ylabel(u'y', fontproperties=prop, rotation=0, labelpad=10)
    ax.set_zlabel(u'z', fontproperties=prop, rotation=0, labelpad=10)

    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    
def plot_S_bloch_regions(ax, barrier=False, init=0):
    state_size = 300
    # Plot |0> and |1>
    ax.scatter(0,0,1.05, label='|0>', marker='^', s=state_size, c='black', zorder=100)
    ax.scatter(0,0,-1.05, label='|1>', marker='v', s=state_size, c='black', zorder=100)
    
    # Plot sphere
    phi = np.linspace(0, np.pi, 20)
    theta = np.linspace(0, 2 * np.pi, 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='w', alpha=0.3, linewidth=0)
    
    # Plot initial
    phi = np.linspace(0, 2*np.pi, 20)
    if init == 0: theta = np.linspace(0, 2*np.arccos(np.sqrt(0.9)), 40)
    elif init == 1: theta = np.linspace(np.pi, 2*np.arccos(np.sqrt(0.1)), 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x,y,z, color='g', alpha=0.7)
    
    # Plot unsafe
    if init == 0: theta = np.linspace(np.pi, 2*np.arccos(np.sqrt(0.89)), 40)
    if init == 1: theta = np.linspace(0, 2*np.arccos(np.sqrt(0.11)), 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    ax.plot_surface(x,y,z, color='r', alpha=0.2)
    
    # Plot barrier
    if barrier:
        X,Y = np.meshgrid(np.arange(-1.1,1.1,step=0.1), np.arange(-1.1,1.1,step=0.1))
        # Check value of f
        f = 0.8 * (-1) ** init
        Z = np.array([[f] * len(X[0])] * len(X))
        xp = X
        yp = Y
        zp = Z
        ax.plot_surface(xp,yp,zp, alpha=0.5)

    ax.set_xlabel(u'x', fontproperties=prop, rotation=0, labelpad=10)
    ax.set_ylabel(u'y', fontproperties=prop, rotation=0, labelpad=10)
    ax.set_zlabel(u'z', fontproperties=prop, rotation=0, labelpad=10)

    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    
def get_states(trange, x_0, ham):
    xs = [np.dot(ham(t), x_0) for t in trange]
    z_vec = [x[0] for x in xs]
    o_vec = [x[1] for x in xs]
    return z_vec, o_vec
    
def state_to_bloch(z_vec, o_vec): 
    mags = [np.absolute(e) for e in z_vec]
    angs = [np.angle(e) for e in z_vec]
    thetas = [2*np.arccos(z) for z in mags]
    phis = []
    for i in range(len(o_vec)):
        phis.append(np.angle(o_vec[i]) - angs[i])
    
    xs = [np.sin(p[0]) * np.cos(p[1]) for p in zip(thetas, phis)]
    ys = [np.sin(p[0]) * np.sin(p[1]) for p in zip(thetas, phis)]
    zs = [np.cos(th) for th in thetas]
    return xs, ys, zs, angs
    
def bloch_transformation(trange, x_0, ham, ax):
    z_vec, o_vec = get_states(trange, x_0, ham)
    
    xs, ys, zs, _ = state_to_bloch(z_vec, o_vec)
    
    ax.scatter(xs[0], ys[0], zs[0], s=100, c='r', zorder=10)
    ax.plot(xs,ys,zs,linestyle='--')
    return ax

def rk4(y0, t0, t, h, f):
    ys = [y0]
    ts = [t0]
    for i in range(int(t/h)):
        k1 = f(ys[i], ts[i])
        k2 = f(ys[i] + h*k1/2, ts[i-1] + h/2)
        k3 = f(ys[i] + h*k2/2, ts[i-1] + h/2)
        k4 = f(ys[i] + h*k3, ts[i-1] + h)
        y = ys[i] + 1/6 * h * (k1 + 2*k2 + 2*k3 + k4)
        ys.append(y)
        t = ts[i] + h
        ts.append(t)
    return np.array(ys), np.array(ts)

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = 0
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )