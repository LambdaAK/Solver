
from cmath import *
from sympy import *

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
        x1 = complex((-self.b-sqrt(d))/(2*self.a))
        x2 = complex((-self.b+sqrt(d))/(2*self.a))
        
        return (x1, x2)
                
                
class HomogenousSecondOrderDifferentialEquation:
    '''
    y'' + a * y' + b * y = 0
    '''
    
    
    def __init__(self, a, b):
        self.a = a
        self.b = b
    
    @property
    def characteristicEquation(self):
        return QuadraticEquation(1, self.a, self.b)
    
    
    def solveCharacteristicEquation(self):
        return self.characteristicEquation.solve()
        
    def solve(self):
        c1, c2 = self.solveCharacteristicEquation()
        
        # check if they're imaginary
        # if they're both real
        if c1.imag == 0 and c2.imag == 0:
            # check if the solutions are lineraly independant
            # if they are lineraly independant
            if c1 != c2:
                return (e ** (c1 * x), e ** (c2 * x))
            # if they are lineraly dependant
            else:
                return (e ** (c1 * x), x * e ** (c2 * x))
            
            
        
        # if they're both imaginary
        '''
        e ^ ax (c1 sin(bx) +  c2 cos(bx))
        
        the complex solutions will be a +- bi
        '''
        
        a = c1.real
        b = abs(c1.imag)

        return (e ** (a * x)) * sin(b * x), (e ** (a * x)) * cos(b * x)
    


ee = HomogenousSecondOrderDifferentialEquation(135135135, 1145987612459876)
     

    