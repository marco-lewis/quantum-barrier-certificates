import itertools

from complex import *
from ComplexVector import *

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.integrate import solve_ivp
import matplotlib as mpl
from matplotlib import animation
from z3 import *

import cmasher as cmr
from IPython.display import HTML
from mpl_toolkits.mplot3d import Axes3D

isqrt2 = 1/np.sqrt(2)
had = np.array([[isqrt2, isqrt2], [isqrt2, -isqrt2]])

trange = np.arange(0, 2*np.pi, 0.01)

X = np.array([[0,1],[1,0]])
Y = np.array([[0,1j],[-1j,0]])
Z = np.array([[1,0],[0,-1]])

prop = mpl.font_manager.FontProperties()
prop.set_file('/usr/share/fonts/truetype/tibetan-machine/TibetanMachineUni.ttf')
prop.set_size(14)

save = lambda file: plt.savefig(file, bbox_inches='tight', pad_inches=0)

def ham(op):
    return lambda t: expm(-1j *op * t)

def had_ham(t):
    return expm(-1j * had * t)

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
        ax.plot_surface(xp,yp,zp, alpha=0.2)

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
    

def print_all_models(s):
    models = []
    check = z3.sat
    while check == z3.sat:
        check = s.check()
        print(check)
        if check == z3.sat:
            model = s.model()
            models.append(model)
            # c = []
            # for var in model: c.append(var() == model[var()])
            # s.add(c)
            # s.push()
            # double_check = s.check()
            # if not(double_check) == z3.sat: raise Exception()
            # s.pop()
            print('Model ' + str(len(models)))
            print(model)
            m = []
            for var in model: m.append(Not(var() == model[var()]))
            s.add(And(m))
    if models == []:
        print(s.unsat_core())
        raise Exception("No models found!")
    return models

def get_real_from_model(model, var, prec=2):
    return round(float(model[var].numerator_as_long())/float(model[var].denominator_as_long()),prec)

def generate_term_powers(k, num_complex_variables):
    x = range(k + 1)
    l = [list(p) for p in itertools.product(x, repeat=2*num_complex_variables)]
    less_than_k = filter(lambda c: sum(c) <= k, l)
    return list(less_than_k)


def round_sympy_expr(expr):
    ret_expr = expr
    for a in sym.preorder_traversal(expr):
        if isinstance(a, sym.Float):
            ret_expr = ret_expr.subs(a, round(a, 2))
    return ret_expr


# Edit of https://stackoverflow.com/a/38980538/19768075
def sympy_to_z3(sympy_var_list, sympy_exp):
    'convert a sympy expression to a z3 expression. This returns (z3_vars, z3_expression)'

    z3_vars = []
    z3_var_map = {}

    for var in sympy_var_list:
        name = var.name
        z3_var = Complex(name)
        z3_var_map[name] = z3_var
        z3_vars.append(z3_var)

    result_exp = _sympy_to_z3_rec(z3_var_map, sympy_exp)

    return z3_vars, result_exp

def _sympy_to_z3_rec(var_map, e):
    'recursive call for sympy_to_z3()'
    core = sym.core

    rv = None

    if not isinstance(e, core.Expr):
        raise RuntimeError("Expected sympy Expr: " + repr(e))

    if isinstance(e, core.Symbol):
        rv = var_map.get(e.name)

        if not(is_complex(rv)):
            raise RuntimeError("No var was corresponds to symbol '" + str(e) + "'")

    elif isinstance(e, core.Number):
        rv = float(e)
    elif isinstance(e, core.Mul):
        rv = _sympy_to_z3_rec(var_map, e.args[0])

        for child in e.args[1:]:
            rv *= _sympy_to_z3_rec(var_map, child)
    elif isinstance(e, core.Add):
        rv = _sympy_to_z3_rec(var_map, e.args[0])

        for child in e.args[1:]:
            rv += _sympy_to_z3_rec(var_map, child)
    elif isinstance(e, core.Pow):
        term = _sympy_to_z3_rec(var_map, e.args[0])
        exponent = _sympy_to_z3_rec(var_map, e.args[1])

        if exponent == 0.5:
            # sqrt
            rv = Sqrt(term)
        else:
            rv = term**exponent
    elif isinstance(e, sym.functions.elementary.complexes.conjugate ):
        rv = _sympy_to_z3_rec(var_map, e.args[0]).conj()
    
    if not(rv):
        raise RuntimeError("Type '" + str(type(e)) + "' is not yet implemented for convertion to a z3 expresion. " + \
                            "Subexpression was '" + str(e) + "'.")

    return rv