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
# Reads HDF5 fields and creates XDMF for visualization
# ---------------------------------------------------------------------

import dolfin as d
import os.path as path

# *** MESH INPUT *** #
meshfile = '../00_meshes/probe_with_mpp_and_porous.xml.gz'

# *** GEOMETRIC PARAMETERS *** #
mm = 0.001  # Milli-meter units

# *** IO PARAMS *** #
inputfile = 'results/%dHz/disp_fields.h5'
frq = 800

# Directory where results are written. End with '/'
results_path = path.abspath('./results/')

# Read HDF5 and write XDMF
mesh = d.Mesh(meshfile)
h5in = d.HDF5File(mesh.mpi_comm(), inputfile % frq, "r")
xdmf_out = d.XDMFFile(path.join(results_path, 'disp_fields3.xdmf'))
xdmf_out.parameters['rewrite_function_mesh'] = False
xdmf_out.parameters['functions_share_mesh'] = True

RTFE = d.FiniteElement("RT", mesh.ufl_cell(), 1)
RT = d.FunctionSpace(mesh, RTFE)
u_re = d.Function(RT)
u_im = d.Function(RT)
u_re.rename("u_re", "0")
u_im.rename("u_im", "0")

h5in.read(u_re, f"/u_sc_re/vector_0")
h5in.read(u_im, f"/u_sc_im/vector_0")

xdmf_out.write(u_re, 0)
xdmf_out.write(u_im, 0)

h5in.close()
xdmf_out.close()