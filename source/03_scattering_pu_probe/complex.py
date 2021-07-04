# ------------------------------------------------------------------ #
#                         ╦═╗╔═╗╔╦╗╔═╗╔═╗╔═╗
#                         ╠╦╝║ ║║║║╚═╗║ ║║  
#                         ╩╚═╚═╝╩ ╩╚═╝╚═╝╚═╝
#  Reduced Order Modelling, Simulation, Optimization of Coupled Systems 
#                             2017-2021
#
#  Authors : 
#  Ashwin Nayak, Andres Prieto, Daniel Fernandez Comesana
#
#  Disclaimer :
#  In downloading this SOFTWARE you are deemed to have read and agreed 
#  to the following terms: This SOFTWARE has been designed with an  
#  exclusive focus on civil applications. It is not to be used for any  
#  illegal, deceptive, misleading or unethical purpose or in any  
#  military applications. This includes ANY APPLICATION WHERE THE USE  
#  OF THE SOFTWARE MAY RESULT IN DEATH, PERSONAL INJURY OR SEVERE  
#  PHYSICAL OR ENVIRONMENTAL DAMAGE. Any redistribution of the software 
#  must retain this disclaimer. BY INSTALLING, COPYING, OR OTHERWISE 
#  USING THE SOFTWARE, YOU AGREE TO THE TERMS ABOVE. IF YOU DO NOT  
#  AGREE TO THESE TERMS, DO NOT INSTALL OR USE THE SOFTWARE.
#
#  Acknowledgements:
#  The ROMSOC project has received funding from the European Union’s 
#  Horizon 2020 research and innovation programme under the Marie 
#  Skłodowska-Curie Grant Agreement No. 765374.
# ------------------------------------------------------------------- #


# -----------------------------------------------------------------
# Definition of a complex class to handle complex-value functions
# -----------------------------------------------------------------
class Complex(object):

    # Init Complex object
    def __init__(self,real_part, imag_part):
        self.real=real_part
        self.imag=imag_part

    # Overload evaluation function
    def __call__(self, x):
        return self.real(x) + 1j * self.imag(x)
    
    # Overload query operator
    def __getitem__(self, key):
        return Complex(self.real[key], self.imag[key])
        
    # Overload addition operator
    def __add__(self, rhs):
        '''Adds complex numbers'''
        return Complex(self.real + rhs.real, self.imag + rhs.imag)
    
    # Overload subtraction operator
    def __sub__(self, rhs):
        '''Subtracts complex numbers'''
        return Complex(self.real - rhs.real, self.imag - rhs.imag)
    
    # Overload Multiplication operator
    def __mul__(self, rhs):
        '''Multiplies complex numbers'''
        return Complex(prod_re(self, rhs), prod_im(self, rhs))
    
    # Overload Division operator
    def __truediv__(self, rhs):
        '''Multiplies complex numbers'''
        return Complex(div_re(self, rhs), div_im(self, rhs))
    

# -------------------------------------------------------------------------------
# Useful functions to handle complex-value functions
# -------------------------------------------------------------------------------

# Define real and imaginary part of the inner-(conjugated) product of two complex-valued functions
def dot_re(f, g):
    return f.real*g.real + f.imag*g.imag
    
def dot_im(f, g):
    return f.imag*g.real - f.real*g.imag 

# Define real and imaginary part of the (unconjugated) product of two complex-valued functions
def prod_re(f, g):
    return f.real*g.real - f.imag*g.imag
    
def prod_im(f, g):
    return f.imag*g.real + f.real*g.imag 

# Define real and imaginary part of the (unconjugated) division of two complex-valued functions
def div_re(f, g):
    return dot_re(f, g)/dot_re(g, g)
    
def div_im(f, g):
    return dot_im(f, g)/dot_re(g, g)