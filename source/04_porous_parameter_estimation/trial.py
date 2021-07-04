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

import h5py
import jax.numpy as np

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class Trial(dict):
    """
    Trial class handling data with the measurements, dict like access
    """
    # Override dict dunder methods
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
    
    # To prevent overload of output
    def __str__(self):
        return 'Trial class object'
    def __repr__(self):
        return 'Trial class object'


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def read_trials(filename=None):
    """
    Reads HDF5 Trial datasets
    """
    if filename is None:
        raise ValueError('Input compatible HDF5file!'
               + "\n Use \x1b[31mpre_process_data.py\x1b[0m to create compatible datasets")
    
    dataset = h5py.File(filename, 'r')
    
    trials = []
    for itrial in range(dataset['ntrials'][()]):
        # Filter data to freq range
        trial = Trial({
            'name' : str(dataset[f'{itrial}']['name'][()]),
            'material' : str(dataset[f'{itrial}']['material'][()]),
            'method' : str(dataset[f'{itrial}']['method'][()]),
            'sample_thickness' : dataset[f'{itrial}']['sample_thickness'][()],
            'probe_distance' : dataset[f'{itrial}']['probe_distance'][()],

            'freqs' : np.array(dataset[f'{itrial}']['freqs']),
            'impedance' : (np.array(dataset[f'{itrial}']['imp_re'])
                          + 1j * np.array(dataset[f'{itrial}']['imp_im'])),
            'reflection' : (np.array(dataset[f'{itrial}']['ref_re'])
                          + 1j * np.array(dataset[f'{itrial}']['ref_im'])),
            'absorption' : np.array(dataset[f'{itrial}']['absorp']),
        })
        trials.append(trial)
    
    dataset.close()
    
    return trials