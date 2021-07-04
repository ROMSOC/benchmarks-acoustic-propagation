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

# Implements the Johnson-Champoux-Allard-Lafarge (JCAL) model to return the porous material properties

# ------------------------------------------------------------------- #


from cmath import sqrt


# Parameters (typical values borrowed from Allard's book)
rho0 = 1.213          # Fluid density
alpha_inf = 1.06      # Static Tortuosity
phi = 0.94            # Porosity
eta = 1.84e-5         # Dynamic Viscosity
sigma = 4e4           # Flow Resistivity
Lmbda = 56e-6         # Viscous Characteristic Length
Lmbda_prime = 110e-6  # Thermal Characteristic Length
Pr = 0.713            # Prandtl Number
k0_prime = eta/sigma  # Static Thermal Permeability
gamma = 1.4           # Ratio of Specific heats
Pi0 = 1.10132e5       # Fluid Equilibrium Pressure


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def rho(w):
    """
    Returns equivalent density of the porous material as a function of frequency
    :param w: Angular Frequency (omega = 2 * pi * frequency)
    :return: Equivalent Density value as a Complex object
    """
    term1 = sqrt(1 + 1j * (4 * alpha_inf**2 * eta * rho0 * w) / (sigma**2 * Lmbda**2 * phi**2))
    rho_eq = (rho0 * alpha_inf / phi) * (1 - 1j * (sigma * phi / (w * rho0 * alpha_inf) * term1))

    return rho_eq


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def bulk_modulus(w):
    """
    Returns equivalent bulk modulus of the porous material as a function of frequency
    :param w: Frequency (omega = 2 * pi * frequency)
    :return: Equivalent Bulk Modulus as a Complex object
    """
    term1 = sqrt(1.0 + 1j * (4 * k0_prime**2 * rho0 * w * Pr) / (eta * Lmbda_prime**2 * phi**2))
    term2 = 1.0 - 1j * (eta * phi / (rho0 * k0_prime * w * Pr)) * term1
    k_eq = (gamma * Pi0 / phi) / (gamma - (gamma - 1) * term2**-1)

    return k_eq
