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
from math import sqrt


def mpp_impedance(w):
    """
    Computes the Impedance of a Micro-Perforated Plate
    in accordance with DY Maa (1987): https://doi.org/10.3397/1.2827694
    :param w: Angular Frequency
    :return: complex-type impedance
    """
    eta = 1.789e-5  # Coefficient of viscosity
    rho = 1.213      # Density of Air
    c = 343.0       # Speed of sound

    p = 0.57        # Ratio of perforated area
    t = 50e-6       # Thickness of MPP
    d = 77e-6       # Diameter of orifice

    x = d * sqrt((w*rho)/(4*eta))

    relative_impedance = (((32*eta*t)/(p*rho*c*d*d)) * (sqrt(1 + x*x/32) + sqrt(x*d/(4*t)))
                          + 1j * (w*t/(p*c)) * (1 + 1/sqrt(9+x*x/2) + 0.85*d/t))
    impedance = relative_impedance * rho * c

    return impedance

#print(mpp_impedance(2*3.141529625*800))
