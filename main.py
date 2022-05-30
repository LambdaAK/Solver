from sympy import *
from cmath import *

x = Symbol('x')

class System:
    '''
    a + b = c
    d + e = f
    '''
    def __init__(self, a, b, c, d, e, f):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        
    def solve(self):
        denominator_matrix_determinant = self.a * self.e - self.b * self.d
        numerator_x_matrix_determinant = self.c * self.e - self.f * self.b
        numerator_y_matrix_determinant = self.a * self.f - self.c * self.d
        
        x = numerator_x_matrix_determinant / denominator_matrix_determinant
        y = numerator_y_matrix_determinant / denominator_matrix_determinant
        
        return {
            'x': x,
            'y': y
        }
        


class LinearDifferentialEquation:
    '''
    dy/dx + py = q
    '''
    def __init__(self, p, q):
        self.p = p
        self.q = q
        
    def solve(self):
        '''
        solve:
        
            mew = e ^ integral(p dx)
            
            d/dx(mew * y) = mew * q
            
            mew * y = integral(mew * q) dx
            
            y = (integral(mew * q) dx) / mew
        '''
        
        mew = e ** integrate(self.p, x)
        
        return (integrate(mew * self.q, x) / mew)
        
        
    
class QuadraticEquation:
    def __init__(self, a, b , c):
        self.a = a
        self.b = b
        self.c = c
        
        
        if a == 0:
            raise Exception('a is zero')
        
        
    def solve(self):
        d = (self.b**2) - (4*self.a*self.c)
        x1 = (-self.b-sqrt(d))/(2*self.a)  
        x2 = (-self.b+sqrt(d))/(2*self.a)
        
        return (x1, x2)
                
                
e = QuadraticEquation(1, 0, 1)

print(e.solve())

