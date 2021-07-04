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

import numpy as np


def compute(inc_mag, w, c, rho_f, rho_p, bulk_mod, r1, r2, r3):
    """
    Computes Coefficients of the pressure field exact solution
    :param g: Boundary condition parameter (constant complex type)
    :param w: Angular Frequency (w = 2 * pi * frequency)
    :param c: Sound Speed
    :param rho_f: Fluid density
    :param rho_p: Porous density
    :param bulk_mod: Porous Bulk Modulus
    :param r1: Radius of Sphere
    :param r2: Inner radius of porous layer
    :param r3: Outer radius of porous layer
    :return: Amplitudes of Exact Pressure field (complex array of size 5)
    """

    # Wave number of fluid and porous layers
    k1 = k3 = w / c
    k2 = np.sqrt(w**2 * rho_p / bulk_mod)

    # Define some functions for ease
    f1 = lambda k, r: np.exp(-1j * k * r) / r
    f2 = lambda k, r: np.exp(1j * k * r) / r
    df1 = lambda k, r: f1(k, r) * (-1.0 - 1j * k * r) / r
    df2 = lambda k, r: f2(k, r) * (-1.0 + 1j * k * r) / r

    # Define 5x5 linear system of equations
    A = np.array([[df1(k1, r1), df2(k1, r1), 0.0, 0.0, 0.0],
                  [f1(k1, r2), f2(k1, r2), -f1(k2, r2), -f2(k2, r2), 0.0],
                  [rho_p * df1(k1, r2), rho_p * df2(k1, r2), -rho_f * df1(k2, r2), -rho_f * df2(k2, r2), 0.0],
                  [0.0, 0.0, f1(k2, r3), f2(k2, r3), -f2(k3, r3)],
                  [0.0, 0.0, rho_f * df1(k2, r3), rho_f * df2(k2, r3), -rho_p * df2(k3, r3)]], dtype=complex)
    b = np.zeros(5, dtype=complex)
    b[0] = - inc_mag * df1(k1, r1)
    coefficients = np.linalg.solve(A, b)

    return coefficients
