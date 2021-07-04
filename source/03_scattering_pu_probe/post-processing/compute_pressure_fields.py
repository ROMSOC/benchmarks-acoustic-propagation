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


# Computes pressure fields from displacement field outputs of PU probe
# ---------------------------------------------------------------------

import dolfin as d
import os.path as path
import complex as C
import numpy as np

# MPP and Porous model modules
# from mpp_impedance import mpp_impedance
import porous_jcal as porous

# Timing modules
from time import localtime, strftime

# Bash script
import sys

# -------------------------------------------------------------------------------------------------
# PARAMETERS
print('[*] Reading Parameters')

# *** MESH INPUT *** #
meshfile = 'probe_with_mpp_and_porous.xml.gz'

# *** GEOMETRIC PARAMETERS *** #
mm = 0.001  # Milli-meter units

# # Probe Params
# Radius of Probe
r_probe = 6.3*mm
# Probe head dimensions
z0 = 28.52*mm
z1 = 44.475*mm

# Porous layer params
t_por = 8.5*mm
h_por = 70*mm

r_0 = r_probe       # Inner radius of porous layer
r_1 = r_0 + t_por   # Outer radius of porous layer
z_0 = -25.53*mm     # Z-coordinate of porous base
z_1 = z_0 + h_por   # Z-coordinate max

# # Domain params
L = [0.02, 0.02, 0.06]
Lpml = 0.015

# *** ACOUSTIC PARAMS *** #
# # Acoustic Fluid Properties
# Sound Speed
c = 343.
# Fluid density
rho = 1.213

# Frequency
freqs = [20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000]

# *** IO PARAMS *** #
inputfile = 'results/%dHz/disp_fields.h5'

# Directory where results are written. End with '/'
results_path = path.abspath('./results/')
# -------------------------------------------------------------------------------------------------
# SOLVER CONFIGURATION

# Mesh sub-domain tolerances
tol_s = d.DOLFIN_EPS_LARGE
tol_l = 1e-7

# Dolfin configuration and optimization
d.parameters["form_compiler"]["optimize"] = True
d.parameters["form_compiler"]["cpp_optimize"] = True
d.parameters['linear_algebra_backend'] = 'PETSc'
# d.set_log_level(d.LogLevel.DEBUG)

# -------------------------------------------------------------------------------------------------
# PROBLEM SETUP
print('[*] Setting up objects acc to parameters')

# Read Mesh
mesh = d.Mesh(meshfile)
mesh.coordinates()[:] *= mm  # Mesh file is in milli-meter units

# *** SUB-DOMAINS *** #
# Fluid domain marker  : 0
# PML domain marker    : 1
# Porous domain marker : 2

sd_marker = {"fluid": 0, "pml": 1, "porous": 2}

# Extract all cell markers
subdomains = d.MeshFunction("size_t", mesh, mesh.topology().dim())
subdomains.rename("Subdomain", "0")

# Expression-based specifications for marking cells in each sub-domain
pml_layer = d.CompiledSubDomain("fabs(x[0])>Lx-tol || fabs(x[1])>Ly-tol || fabs(x[2])>Lz-tol",
                                Lx=L[0], Ly=L[1], Lz=L[2], tol=tol_s)
porous_layer = d.CompiledSubDomain(""" x[0]*x[0] + x[1]*x[1] > r0*r0 - tol 
                                       &&  x[0]*x[0] + x[1]*x[1] < r1*r1 + tol
                                       && x[2] > z0-tol && x[2] < z1+tol""",
                                   r0=r_0, r1=r_1, z0=z_0, z1=z_1, tol=tol_l)

# Mark cells in sub-domains
subdomains.set_all(sd_marker["fluid"])
pml_layer.mark(subdomains, sd_marker["pml"])
porous_layer.mark(subdomains, sd_marker["porous"])

# *** NORMALS *** #
# Extract face normal vectors
# NOTE: All normals point outwards of the enclosed domain
n_vec = d.FacetNormal(mesh)

# *** FUNCTION SPACE *** #
# Define function space of first-order Raviart-Thomas element
RTFE = d.FiniteElement("RT", mesh.ufl_cell(), 1)
RT = d.FunctionSpace(mesh, RTFE)

# Define function space for pressure field
DGFE = d.FiniteElement("DG", mesh.ufl_cell(), 0)
DG0 = d.FunctionSpace(mesh, DGFE)
DG0Complex = d.FunctionSpace(mesh, DGFE*DGFE)

# *** I/O *** #
# Field outputs in HDF5
outfile = path.join(results_path, 'pres_fields.h5')
h5out = d.HDF5File(mesh.mpi_comm(), outfile, "w")

# -------------------------------------------------------------------------------------------------
# SOLVER SETUP
dx = d.Measure('dx', domain=mesh, subdomain_data=subdomains)
A_mat = d.PETScMatrix()
b_vec = d.PETScVector()

# -------------------------------------------------------------------------------------------------
# LOOP OVER FREQUENCY
for frq_count, frq in enumerate(freqs):
    print(f'[*] Computing Pressure field for {frq}Hz...  ({frq_count+1}/{len(freqs)})')
    print('[ ][*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))

    # Open input HDF5 file
    h5in = d.HDF5File(mesh.mpi_comm(), inputfile % frq, "r")

    # Read disp field from file
    # Read HDF5 File for the relevant disp fields
    print('[ ][*] Reading Disp fields')
    u_re = d.Function(RT)
    u_im = d.Function(RT)
    h5in.read(u_re, f"/u_tot_re/vector_0")
    h5in.read(u_im, f"/u_tot_im/vector_0")

    # Angular frequency
    omega = 2 * d.DOLFIN_PI * frq

    # Compute Dyn Bulk-modulus
    BulkMod = porous.bulk_modulus(omega)
    CBulkMod = C.Complex(float(BulkMod.real), float(BulkMod.imag))

    # Real part of Pressure
    print('[ ][*] Computing Pressure field')

    # Test and Trial functions
    (w_re, w_im) = d.TestFunction(DG0Complex)
    (Pv_re, Pv_im) = d.TrialFunction(DG0Complex)

    # # Project displacement solution on the DG0 space
    div_u = C.Complex(d.div(u_re), d.div(u_im))

    a = (d.inner(w_re, Pv_re) + d.inner(w_im, Pv_im)) * dx
    rhs = ((d.inner(w_re, - rho * c * c * d.div(u_re))
            + d.inner(w_im, - rho * c * c * d.div(u_im))) * dx(sd_marker["fluid"])
           + (d.inner(w_re, - rho * c * c * d.div(u_re))
              + d.inner(w_im, - rho * c * c * d.div(u_im))) * dx(sd_marker["pml"])
           + (d.inner(w_re, - C.prod_re(CBulkMod, div_u))
              + d.inner(w_im, - C.prod_im(CBulkMod, div_u))) * dx(sd_marker["porous"]))

    A = d.assemble(a, tensor=A_mat)
    b = d.assemble(rhs, tensor=b_vec)

    p = d.Function(DG0Complex)
    solver = d.LUSolver(A, 'mumps')
    solver.solve(p.vector(), b)

    (p_re, p_im) = p.split(True)

    # Output to file
    print('[ ][*] Writing pres field to file.')
    p_re.rename("p_re", "Pressure real")
    p_im.rename("p_im", "Pressure imag")
    h5out.write(p_re, "p_re", frq_count)
    h5out.write(p_im, "p_im", frq_count)

    h5in.close()

h5out.close()
