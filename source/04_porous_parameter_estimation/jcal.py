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
import jax.numpy as np


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class PorousModel:
    """
    Porous Model super class which provides dict like access to params
    """
    def set_params(self, p):
        if isinstance(p, dict):
            for k, v in p.items():
                if k in self.__dict__:
                    self.__dict__[k] = v
                else:
                    raise KeyError(f"Parameter '{k}'  unrecognized!")
        else:
            raise ValueError('Input dict argument')
    
    def get_fitting_params(self):
        raise NotImplementedError

    def set_fitting_params(self, params):
        raise NotImplementedError

    def density(self, w):
        raise NotImplementedError
    
    def bulk_modulus(self, w):
        raise NotImplementedError
    
    def wavenumber(self, w):
        raise NotImplementedError
    
    def characteristic_impedance(self, w):
        raise NotImplementedError
    

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class JCALModel(PorousModel):
    """
    JCAL Model with classic 6 parameters:
    alpha_infty, sigma, phi, lmbda, lmbda_prime, k0_prime
    """
    
    def __init__(self):
        # Porous Parameters
        self.alpha_infty = 1.06    # Static Tortuosity
        self.phi = 0.94            # Porosity
        self.sigma = 15.2e3          # Flow Resistivity
        self.lmbda =  56e-6          # Viscous Characteristic Length
        self.lmbda_prime = 110e-6    # Thermal Characteristic Length
        self.k0_prime = 2.76e-10       # Static Thermal Permeability
        
        # Fluid Parameters
        self.rho_f = 1.213         # Fluid density
        self.c_f = 343.0           # Speed of sound
        self.Pr = 0.713            # Prandtl Number
        self.eta = 1.84e-5         # Dynamic Viscosity or sigma * k0_prime
        self.gamma = 1.4           # Ratio of Specific heats
        self.Pi0 = 1.0132e5       # Fluid Equilibrium Pressure
    
    def get_fitting_params(self):
        return np.array([
            self.alpha_infty, 
            self.phi,
            self.sigma,
            self.lmbda,
            self.lmbda_prime,
            self.k0_prime
            ])
    
    def set_fitting_params(self, params):
        ( 
            self.alpha_infty, 
            self.phi,
            self.sigma,
            self.lmbda,
            self.lmbda_prime,
            self.k0_prime
        ) = params
        pass
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def density(self, w):
        """
        Returns equivalent density of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Density value as a Complex object
        """
        temp = np.sqrt( 1 + 1j * (4 * self.alpha_infty**2 * self.eta * self.rho_f * w) 
                        / (self.sigma**2 * self.lmbda**2 * self.phi**2))
        return ((self.rho_f * self.alpha_infty / self.phi) 
                    * (1 - 1j * self.sigma * self.phi * temp / (w * self.rho_f * self.alpha_infty) ))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def bulk_modulus(self, w):
        """
        Returns equivalent bulk modulus of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Bulk Modulus as a Complex object
        """
        temp1 = np.sqrt( 1.0 + 1j * (4 * self.k0_prime**2 * self.rho_f * w * self.Pr) 
                            / (self.eta * self.lmbda_prime**2 * self.phi**2))
        temp2 = 1.0 - 1j * (self.eta * self.phi 
                        / (self.rho_f * self.k0_prime * w * self.Pr)) * temp1
        return ((self.gamma * self.Pi0 / self.phi) 
                / (self.gamma - (self.gamma - 1) * temp2**-1))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def wavenumber(self, w):
        """
        Returns porous wavenumber 
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        """
        return np.sqrt(self.density(w) * w * w / self.bulk_modulus(w))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def characteristic_impedance(self, w):
        return np.sqrt(self.density(w) * self.bulk_modulus(w) * self.phi * self.phi)
    

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class JCALModel2(PorousModel):
    """
    JCAL Model with reduced params : 
    alpha_infty, sigma*phi, lmbda, lmbda_prime, phi / k0_prime 
    """
    
    def __init__(self):
        # Porous Parameters
        self.alpha_infty = 1.06    # Static Tortuosity
        self.sigma_phi = 15.2e3 * 0.94
        self.lmbda =  56e-6          # Viscous Characteristic Length
        self.lmbda_prime = 110e-6    # Thermal Characteristic Length
        self.phi_div_k0_prime = 0.94 / 2.76e-4
        
        # Fluid Parameters
        self.rho_f = 1.213         # Fluid density
        self.c_f = 343.0           # Speed of sound
        self.Pr = 0.713            # Prandtl Number
        self.eta = 1.84e-5         # Dynamic Viscosity or sigma * k0_prime
        self.gamma = 1.4           # Ratio of Specific heats
        self.Pi0 = 1.0132e5       # Fluid Equilibrium Pressure
    
    def get_fitting_params(self):
        return np.array([
            self.alpha_infty, 
            self.sigma_phi,
            self.lmbda,
            self.lmbda_prime,
            self.phi_div_k0_prime
            ])
    
    def set_fitting_params(self, params):
        ( 
            self.alpha_infty, 
            self.sigma_phi,
            self.lmbda,
            self.lmbda_prime,
            self.phi_div_k0_prime
        ) = params
        pass
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def phi_density(self, w):
        """
        Returns equivalent density of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Density value as a Complex object
        """
        temp = np.sqrt( 1 + 1j * (4 * self.alpha_infty**2 * self.eta * self.rho_f * w) 
                        / (self.sigma_phi**2 * self.lmbda**2))
        return ((self.rho_f * self.alpha_infty ) 
                    * (1 - 1j * self.sigma_phi * temp / (w * self.rho_f * self.alpha_infty) ))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def phi_bulk_modulus(self, w):
        """
        Returns equivalent bulk modulus of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Bulk Modulus as a Complex object
        """
        temp1 = np.sqrt( 1.0 + 1j * (4  * self.rho_f * w * self.Pr) 
                            / (self.eta * self.lmbda_prime**2 * self.phi_div_k0_prime**2))
        temp2 = 1.0 - 1j * (self.eta * self.phi_div_k0_prime
                        / (self.rho_f *  w * self.Pr)) * temp1
        return ((self.gamma * self.Pi0 ) 
                / (self.gamma - (self.gamma - 1) * temp2**-1))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def wavenumber(self, w):
        """
        Returns porous wavenumber 
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        """
        return np.sqrt(self.phi_density(w) * w * w / self.phi_bulk_modulus(w))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def characteristic_impedance(self, w):
        return np.sqrt(self.phi_density(w) * self.phi_bulk_modulus(w))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class JCALModel3(PorousModel):
    """
    JCAL Model with reduced and SCALED params : 
    alpha_infty, sigma*phi, lmbda, lmbda_prime, phi / k0_prime 
    """
    
    def __init__(self):
        # Porous Parameters
        self.alpha_infty = 1.06    # Static Tortuosity
        self.sigma_phi = 14.0       # x 10**3
        self.lmbda =  56.0          # x 10**-6 Viscous Characteristic Length
        self.lmbda_prime = 110.0    # x 10**-6 Thermal Characteristic Length
        self.phi_div_k0_prime = 1.0 # x 10**10 Static Thermal Permeability
        
        # Fluid Parameters
        self.rho_f = 1.213         # Fluid density
        self.c_f = 343.0           # Speed of sound
        self.Pr = 0.713            # Prandtl Number
        self.eta = 1.84e-5         # Dynamic Viscosity or sigma * k0_prime
        self.gamma = 1.4           # Ratio of Specific heats
        self.Pi0 = 1.0132e5       # Fluid Equilibrium Pressure
    
    def get_fitting_params(self):
        return np.array([
            self.alpha_infty, 
            self.sigma_phi,
            self.lmbda,
            self.lmbda_prime,
            self.phi_div_k0_prime
            ])
    
    def set_fitting_params(self, params):
        (   self.alpha_infty, 
            self.sigma_phi,
            self.lmbda,
            self.lmbda_prime,
            self.phi_div_k0_prime
        ) = params
        pass
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def phi_density(self, w):
        """
        Returns equivalent density of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Density value as a Complex object
        """
        temp = np.sqrt( 1 + 1j * (4 * self.alpha_infty**2 * self.eta * self.rho_f * w) 
                        / (self.sigma_phi**2  * self.lmbda**2 * 1e-6))
        return ((self.rho_f * self.alpha_infty ) 
                    * (1 - 1j * self.sigma_phi * 1e3 * temp / (w * self.rho_f * self.alpha_infty) ))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def phi_bulk_modulus(self, w):
        """
        Returns equivalent bulk modulus of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Bulk Modulus as a Complex object
        """
        temp1 = np.sqrt( 1.0 + 1j * (4  * self.rho_f * w * self.Pr) 
                            / (self.eta * self.lmbda_prime**2 * 1e8 * self.phi_div_k0_prime**2 ))
        temp2 = 1.0 - 1j * (self.eta * self.phi_div_k0_prime * 1e10
                        / (self.rho_f *  w * self.Pr)) * temp1
        return ((self.gamma * self.Pi0 )
                / (self.gamma - (self.gamma - 1) * temp2**-1))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def wavenumber(self, w):
        """
        Returns porous wavenumber 
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        """
        return np.sqrt(self.phi_density(w) * w * w / self.phi_bulk_modulus(w))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def characteristic_impedance(self, w):
        return np.sqrt(self.phi_density(w) * self.phi_bulk_modulus(w))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class JCALModel4(PorousModel):
    """
    JCAL Model with reduced params (3 param Heroshenkov) : 
    Porosity (phi), Median pore size (s), Std deviation of Pore size (sigma_s)
    """
    
    def __init__(self):
        # Porous Parameters
        self.alpha_infty = 1.06    # Static Tortuosity
        self.phi = 0.94            # Porosity
        self.sigma = 15.2e3          # Flow Resistivity
        self.lmbda =  56e-6          # Viscous Characteristic Length
        self.lmbda_prime = 110e-6    # Thermal Characteristic Length
        self.k0_prime = 2.76e-10       # Static Thermal Permeability

        self.s = 1e-6 # Median Pore Size
        self.sigma_s = 2e-7 # Std deviation in Pore Size
        
        # Fluid Parameters
        self.rho_f = 1.213         # Fluid density
        self.c_f = 343.0           # Speed of sound
        self.Pr = 0.713            # Prandtl Number
        self.eta = 1.84e-5         # Dynamic Viscosity or sigma * k0_prime
        self.gamma = 1.4           # Ratio of Specific heats
        self.Pi0 = 1.0132e5       # Fluid Equilibrium Pressure
    
    def get_fitting_params(self):
        return np.array([ 
            self.phi,
            self.s,
            self.sigma_s,
            ])
    
    def set_fitting_params(self, params):
        (   self.phi,
            self.s,
            self.sigma_s,
        ) = params
        self.alpha_infty = self.func_alpha_infty()
        self.sigma = self.func_sigma()
        self.lmbda = self.func_lmbda()
        self.lmbda_prime = self.func_lmbda_prime()
        self.k0_prime = self.func_k0_prime()
        pass
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def func_alpha_infty(self):
        return np.exp(4*(self.sigma_s*np.log10(2))**2)
    
    def func_sigma(self):
        return 8*self.eta*self.s**2 * np.exp(10*(self.sigma_s*np.log10(2))**2) / self.phi
    
    def func_lmbda(self):
        return self.s * np.exp(-2.5*(self.sigma_s*np.log10(2))**2)

    def func_lmbda_prime(self):
        return self.s * np.exp(1.5*(self.sigma_s*np.log10(2))**2)
    
    def func_k0_prime(self):
        return self.s**2*self.phi*np.exp(2*(self.sigma_s*np.log10(2))**2)/8

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def density(self, w):
        """
        Returns equivalent density of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Density value as a Complex object
        """
        temp = np.sqrt( 1 + 1j * (4 * self.alpha_infty**2 * self.eta * self.rho_f * w) 
                        / (self.sigma**2 * self.lmbda**2 * self.phi**2))
        return ((self.rho_f * self.alpha_infty / self.phi) 
                    * (1 - 1j * self.sigma * self.phi * temp / (w * self.rho_f * self.alpha_infty) ))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def bulk_modulus(self, w):
        """
        Returns equivalent bulk modulus of the porous material as a function of frequency
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        :return: Equivalent Bulk Modulus as a Complex object
        """
        temp1 = np.sqrt( 1.0 + 1j * (4 * self.k0_prime**2 * self.rho_f * w * self.Pr) 
                            / (self.eta * self.lmbda_prime**2 * self.phi**2))
        temp2 = 1.0 - 1j * (self.eta * self.phi 
                        / (self.rho_f * self.k0_prime * w * self.Pr)) * temp1
        return ((self.gamma * self.Pi0 / self.phi) 
                / (self.gamma - (self.gamma - 1) * temp2**-1))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def wavenumber(self, w):
        """
        Returns porous wavenumber 
        :param w: Angular Frequency (omega = 2 * pi * frequency)
        """
        return np.sqrt(self.density(w) * w * w / self.bulk_modulus(w))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def characteristic_impedance(self, w):
        return np.sqrt(self.density(w) * self.bulk_modulus(w) * self.phi * self.phi)