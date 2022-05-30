
from cmath import *
from sympy import *

x = Symbol('x')
c1 = Symbol('c1')
c2 = Symbol('c2')




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
        
        the complex solutions will be in the form a +- bi
        '''
        
        a = c1.real
        b = abs(c1.imag)

        return (e ** (a * x)) * sin(b * x), (e ** (a * x)) * cos(b * x)
    

class NonhomogeneousSecondOrderDifferentialEquation(HomogenousSecondOrderDifferentialEquation):
    '''
    Uses variation of parameters to solve for the particular solution to the differential equation and
    Uses a characteristic equation and the quadratic formula to solve for the characteristic solutions to the differential equation
    '''
    
    def __init__(self, a, b, r):
        super().__init__(a, b)
        self.r = r
    
    
    
    
    def findParticularSolution(self):
        '''
        u'y1 + v'y2 = 0
        u'y1' + v'y2' = r
        
        
        The particular solution is uy1 + vy2 where 
            u and v are the solutions to the system of equations defined above
            y1 and y2 are the general characteristic solutions to the corresponding homogenous differential equation where r = 0
        '''
        
        y1, y2 = super().solve()
        
        # create the system of equations
        '''
        u_prime = Symbol('u\'')
        v_prime = Symbol('v\'')
        '''
        y1_prime = diff(y1, x)
        y2_prime = diff(y2, x)
        
     
        
        
      
        
        # cramer's rule
        '''
        demoninator
        
        y1       y2
        y1_prime y2_prime
        
        
        numerator u'
        
        
        0       y2
        r       y2_prime
        
        numerator v'
        
        y1       0
        y1_prime r
        '''
        
        
        
        denominator = Matrix([[y1, y2], [y1_prime, y2_prime]])
        numerator_u_prime = Matrix([[0, y2], [self.r, y2_prime]])
        numerator_v_prime = Matrix([[y1, 0], [y1_prime, self.r]])
        
        
        u_prime = numerator_u_prime.det() / denominator.det()
        v_prime = numerator_v_prime.det() / denominator.det()
        
        
        u = integrate(u_prime, x)
        v = integrate(v_prime, x)
        
        
        return u * y1 + v * y2

    def solve(self):
        y1, y2 = super().solve()
        
        yp = self.findParticularSolution()
        
        return c1 * y1 + c2 * y2 + yp
        
        
        
        
        
        
        
        
    
test = NonhomogeneousSecondOrderDifferentialEquation(-2, 1, x)

print(test.solve())