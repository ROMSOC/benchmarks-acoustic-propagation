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
# Computes total displacement field at a specific coordinate
#  - of PU probe upon incidence by a plane wave

import dolfin as d
import json

# PARAMETERS
print('[*] Reading Parameters')

# Mesh 
meshfile = 'probe_with_mpp_and_porous.xml.gz'
mm = 0.001  # Mesh units

# Frequency
freqs_list = [20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000]

# Input File with displacement fields
inputfile = 'results/%dHz/disp_fields.h5'

# Evaluation Point coordinates
p1 = [0.0, 0.0, 0.035]
points_list = [p1]

# CONFIGURATION OPTIONS
# Fields to read from input file
# TODO: Need to fix this by iterating over fields
finput_field = 'u_tot'

# Output file name
outfile = "point_evaluations_disp.json"
# Dolfin configuration options
#d.set_log_level(d.LogLevel.DEBUG)

# -------------------------------------------------------------------------------------------------
# SETUP
print('[*] Setting up objects acc to parameters')

# Read Mesh
mesh = d.Mesh(meshfile)
mesh.coordinates()[:] *= mm  # Mesh file is in milli-meter units

# *** FUNCTION SPACE *** #
# Define function space of first-order Raviart-Thomas element
RTFE = d.FiniteElement("RT", mesh.ufl_cell(), 1)
RT = d.FunctionSpace(mesh, RTFE)

# Init results data structure
points_dict = {f'point{point_count:d}': {'coordinate': point} for point_count, point in enumerate(points_list)}

# Loop over frequencies and read specific files
for frq_count, frq in enumerate(freqs_list):
    print(f'[*] Reading data for {frq}Hz...  ({frq_count+1}/{len(freqs_list)})')

    # Open input HDF5 file
    h5id = d.HDF5File(mesh.mpi_comm(), inputfile % frq, "r")

    # Read disp field from file
    field_re = d.Function(RT)
    field_im = d.Function(RT)
    h5id.read(field_re, f"/{finput_field}_re/vector_0")
    h5id.read(field_im, f"/{finput_field}_im/vector_0")

    # Iterate over points
    for point_count, point in enumerate(points_list):
        # Evaluate function at point
        field_re_pt = field_re(point)
        field_im_pt = field_im(point)

        # Write to result data-structure
        points_dict[f"point{point_count:d}"][f"freq{frq_count:d}"] = {"value": frq,
                                                                      f"{finput_field}_re": field_re_pt.tolist(),
                                                                      f"{finput_field}_im": field_im_pt.tolist()}

    h5id.close()

# Write result to file
print(f'Results written to {outfile}.')
with open(outfile, "w") as file:
    json.dump(points_dict, file, indent=4, sort_keys=True)
