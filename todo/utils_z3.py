from complex import *
from ComplexVector import *
import sympy as sym
from z3 import *

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