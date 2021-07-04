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
import jcal as jc
import trial as tr

import jax
import jax.numpy as np
import scipy as sp
import json

from time import localtime, strftime

# Config
jax.config.update('jax_platform_name', 'cpu')
jax.config.update('jax_enable_x64', True)
np.set_printoptions(precision=3)

# Data
filename = 'data/experiments.data' 

# Initial Guess for params
p0 = np.array(
        [1.2, 12.0, 100.0, 100.0, 0.1], 
        dtype=np.float64) 


# Bounds
bounds = [(1.1, 1.3), (1, 50), (1, 400), (50, 1000), (0.01, 1.0)]

# Read data
trials = tr.read_trials(filename)
print('Number of datasets : ', len(trials))

# JCAL Model
jcal = jc.JCALModel3()
jcal.rho_f = 1.2047
jcal.c_f = 343.149

# Impedance Error 
def error_imp(p):
    """
    Overall Impedance Error (w.r.to log(freq) scale)
    """
    jcal.set_fitting_params(p)
    error = 0.
    for trial in trials:
        omega_vals = 2*np.pi*trial['freqs']
        data = trial['impedance']*jcal.rho_f*jcal.c_f
        s_imp = (jcal.characteristic_impedance(omega_vals)
                 / np.tanh(1j*jcal.wavenumber(omega_vals)*trial['sample_thickness']))
        error = error + np.linalg.norm((data - s_imp)/trial['freqs']) / np.linalg.norm(data/trial['freqs'])
    return error/len(trials)

# Cost Function 
cost_func = jax.jit(error_imp, backend='cpu')
jacobian = jax.grad(cost_func)

# Minimization
Minimas= {}

print('# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print('[*] Minimization with TNC Method and bounds')
print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))

Minima1 = sp.optimize.minimize(
                cost_func, 
                p0, 
                method='TNC',
                bounds=bounds,
                options={
                    'eta' : 1,
                    'maxiter':75000,
                    'ftol': 1e-20,
                    # 'disp': True 
                    }
                )

print(f"""
        Cost Function   : {Minima1.fun}
        N Iters         : {Minima1.nit}
        N Fun Evals     : {Minima1.nfev}
        Init Params     : {p0}
        Solution Params : {Minima1.x} 
       """)
Minimas['TNC'] = (Minima1.x).tolist()


print('# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print('[*] Minimization with SLSQP Method with provided Jacobian and Bounds')
print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))

Minima2 = sp.optimize.minimize(
                cost_func, 
                p0, 
                method='SLSQP',
                jac=jacobian,
                # hess=hessian, 
                bounds=bounds,
                options={
                    'maxiter':75000,
                    'ftol' : 1e-20
                    },
                tol=1e-20
                )

print(f"""
        Cost Function   : {Minima2.fun}
        N Iters         : {Minima2.nit}
        N Fun Evals     : {Minima2.nfev}
        Init Params     : {p0}
        Solution Params : {Minima2.x} 
       """)

print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))
Minimas['SLSQP2'] = (Minima2.x).tolist()


print('# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
# Scipy Least Squares

# Error Vector Function 
error_vec_size = 0
for trial in trials:
    error_vec_size += trial['freqs'].size 
error_vec = np.empty((error_vec_size,), dtype=np.complex128)

def error_imp_vec(p):
    jcal.set_fitting_params(p)
    counter_start = 0
    norm_vec = error_vec
    for trial in trials:
        omega_vals = 2*np.pi*trial['freqs']
        
        s_imp = (jcal.characteristic_impedance(omega_vals)
                 / np.tanh(1j*jcal.wavenumber(omega_vals)*trial['sample_thickness']))
        counter_end = counter_start + trial['freqs'].size
        norm_vec = jax.ops.index_update(
            norm_vec,
            jax.ops.index[counter_start:counter_end], 
            (trial['impedance']*jcal.rho_f*jcal.c_f - s_imp)/omega_vals, 
            indices_are_sorted= True,
            unique_indices = False
        )
        counter_start = counter_end
    
    return np.abs(norm_vec)


# Cost Function
cost_func2 = jax.jit(error_imp_vec, backend='cpu')
jacobian2 = jax.jacobian(cost_func2)
bounds2 = ([a[0] for a in bounds], [a[1] for a in bounds])


print('[*] Minimization with Trust-Region Least Squares and Bounds')
print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))

Minima3 = sp.optimize.least_squares(
                cost_func2,
                p0,
                jac=jacobian2,
                bounds=bounds2,
                method='trf', # Trust-region Reflective
                ftol=1e-15,
                xtol=1e-15,
                gtol=1e-15,
                x_scale='jac',
                loss='linear',
                max_nfev=50000,
                tr_solver='exact',
                verbose=1,
    )

print(f"""
        Cost Function   : {Minima3.cost}
        N Fun Evals     : {Minima3.nfev}
        Init Params     : {p0}
        Solution Params : {Minima3.x} 
       """)

print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))
Minimas['Least_Squares'] = (Minima3.x).tolist()

print('# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

# Write Minimas to file
with open("fit_params.json", 'w') as file:
    json.dump(Minimas, file, indent=2, sort_keys=True)