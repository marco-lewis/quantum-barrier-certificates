from z3 import *

class InitRegion():
    def __init__(self, data_gen):
        self.data_gen = data_gen
    
    def add_conditions(self, solver,barrier, n=10):
        points = self.get_points(n)
        for x in points:
            bx = barrier.apply(x)
            solver.add(simplify(bx.r <= 0), simplify(bx.i == 0))
    
    def get_points(self, n):
        return [self.data_gen() for i in range(n)]
    
class UnsafeRegion():
    def __init__(self, data_gen):
        self.data_gen = data_gen
    
    def add_conditions(self, solver, barrier, n=10):
        points = self.get_points(n)
        for x in points:
            bx = barrier.apply(x)
            solver.add(simplify(bx.r > 0), simplify(bx.i == 0))
    
    def get_points(self, n):
        return [self.data_gen() for i in range(n)]