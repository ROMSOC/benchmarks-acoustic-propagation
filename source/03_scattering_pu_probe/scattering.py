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
#
# Computes scattered displacement field of PU probe upon incidence by a plane wave
# --------------------------------------------------------------------------------
# - Includes MPP model to account for steel mesh around the probe head.
# - Includes PML model to account for replicating free-field conditions
# - Includes Porous model
#

import dolfin as d
import complex as C
import os.path as path

# MPP / Porous model modules
from mpp_impedance import mpp_impedance
import porous_jcal as porous

# Timing modules
from time import localtime, strftime

# Bash script
import sys

# -------------------------------------------------------------------------------------------------
# PARAMETERS
print('[*] Reading Parameters')

# *** MESH INPUT *** #
meshfile = '../00_meshes/probe_with_mpp_and_porous.xml.gz'

# *** GEOMETRIC PARAMETERS *** #
mm = 0.001  # Milli-meter units in mesh

# # Probe Params
# Radius of Probe
r_probe = 6.3*mm
# Probe head dimensions
z0 = 28.52*mm
z1 = 44.475*mm

# Porous layer params
t_por = 8.5*mm      # Thickness
r_0 = r_probe       # Inner radius of porous layer
r_1 = r_0 + t_por   # Outer radius of porous layer

h_por = 70*mm
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

# # Incident plane wave params
# Amplitude
inc_mag = 1.
# Direction
inc_dir = [1., 0., 0.]  # Direction

# Frequency
freqs = [float(sys.argv[i]) for i in range(1, len(sys.argv))]
freqs_list = [20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000]
freqs_dict = {freqs_list[i]: i for i in range(len(freqs_list))}

# *** IO PARAMS *** #
# Directory where results are written. End with '/'
results_path = path.abspath('./results/')

# -------------------------------------------------------------------------------------------------
# SOLVER CONFIGURATION

# Mesh sub-domain tolerances
tol_s = d.DOLFIN_EPS_LARGE
tol_l = 1e-7

# Diagnostic outputs
loutdiag = False
loutmarkers = False

# IO Configuration
h5_write_mode = "w"

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

# *** BOUNDARIES *** #
# Generic faces         : 0
# Rigid object boundary : 1
# MPP interface         : 2
# PML boundary          : 3
bd_marker = {"generic": 0, "object": 1, "mpp": 2, "pml": 3}

# Extract all face markers
boundaries = d.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.rename("Boundary", "0")

# Expression-based specifications for marking faces
pml_bnd = d.CompiledSubDomain("""on_boundary 
                                 && ( near(fabs(x[0]), Lx+Lpml) 
                                      || near(fabs(x[1]),Ly+Lpml)
                                      || near(fabs(x[2]),Lz+Lpml) 
                                 )""",
                              Lx=L[0], Ly=L[1], Lz=L[2], Lpml=Lpml, tol=tol_s)
object_bnd = d.CompiledSubDomain("on_boundary && fabs(x[0])<Lx+tol && fabs(x[1])<Ly+tol && fabs(x[2])<Lz+tol",
                                 Lx=L[0], Ly=L[1], Lz=L[2], tol=tol_s)
mpp_bnd = d.CompiledSubDomain("""x[2] > z0 && x[2]< z1 
                                 && near(x[0]*x[0] + x[1]*x[1], radius*radius, tol)""",
                              z0=z0, z1=z1, radius=r_probe, tol=tol_l)

# # Mark boundaries
boundaries.set_all(bd_marker["generic"])
pml_bnd.mark(boundaries, bd_marker["pml"])
mpp_bnd.mark(boundaries, bd_marker["mpp"])
object_bnd.mark(boundaries, bd_marker["object"])

# *** NORMALS *** #
# Extract face normal vectors
# NOTE: All normals point outwards of the enclosed domain
n_vec = d.FacetNormal(mesh)

# *** FUNCTION SPACE *** #
# Define function space of first-order Raviart-Thomas element
RTFE = d.FiniteElement("RT", mesh.ufl_cell(), 1)
RT = d.FunctionSpace(mesh, RTFE)

# Define 2D real space to substitute complex space
# NOTE: sub(0) and sub(1) attributes are the subspaces
RTComplex = d.FunctionSpace(mesh, RTFE * RTFE)

# *** I/O *** #
# Field outputs in HDF5
outfile = path.join(results_path, str(int(freqs[0])) + 'Hz', 'disp_fields.h5')
h5out = d.HDF5File(mesh.mpi_comm(), outfile, h5_write_mode)

# Diagnostic Output
if loutdiag:
    # Problem Configuration log
    finfo = open(path.join(results_path, 'INFO.log'), 'w')
    finfo.write('MESH INFO:\n')
    finfo.write('Mesh Name = '+meshfile)
    finfo.write('\nNumber of vertices in mesh = ' + str(mesh.coordinates().shape[0])+'\n')
    finfo.write('h_min = ' + str(mesh.hmin()) + '\n')
    finfo.write('h_max = ' + str(mesh.hmax()) + '\n')
    finfo.write('\nPHYSICAL PARAMETERS : \n')
    finfo.write('Lx, Ly, Lz = '+str(L)+'\nLpml = '+str(Lpml)+'\n')
    finfo.write('\nFLUID PARAMETERS:\n')
    finfo.write('c = '+str(c)+'\nrho = '+str(rho)+'\n')
    finfo.write('\nPROBLEM CONFIG:\n')
    finfo.write('Freqs = '+str(freqs)+'\n')
    finfo.close()


if loutmarkers:
    # Boundary and domain markers to file
    fboundaries = d.File(path.join(results_path, 'BoundaryMarkers.pvd'))
    fboundaries << boundaries
    fdomains = d.File(path.join(results_path, 'DomainMarkers.pvd'))
    fdomains << subdomains


# -------------------------------------------------------------------------------------------------
# SOLVER SETUP
print('[*] Setting up solver objects and expressions')

# *** Variational form objects and expressions ***
# Measures for integrals
# Obtain volume integral measures for each sub-domain
dx = d.Measure('dx', domain=mesh, subdomain_data=subdomains)

# Obtain surface integral measures for boundaries
ds = d.Measure('ds', domain=mesh, subdomain_data=boundaries)

# Surface integral measure for MPP internal face
dS = d.Measure('dS', domain=mesh, subdomain_data=boundaries, subdomain_id=bd_marker["mpp"])

# Zero Expressions
zero = d.Constant("0.")
zero_vec = d.Constant((0., 0., 0.))

# Boundary Condition expressions
# g = u_sc :dot: n = - u_inc :dot: n
# where u_inc : Incident plane wave disp field
cpp_g_re = ("-A*k*k0_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k1_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k2_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_re = d.Expression(cpp_code=cpp_g_re,
                    A=0.,                     # Will be overwritten later
                    k=0.,                     # Will be overwritten later
                    k0_=inc_dir[0], k1_=inc_dir[1], k2_=inc_dir[2],
                    degree=4)

cpp_g_im = ("A*k*k0_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k1_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k2_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_im = d.Expression(cpp_code=cpp_g_im,
                    A=0.,                   # Will be overwritten later
                    k=0.,                   # Will be overwritten later
                    k0_=inc_dir[0], k1_=inc_dir[1], k2_=inc_dir[2],
                    degree=4)

# # RHS
cpp_uinc_re = ("-A*k*k0_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
               "-A*k*k1_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
               "-A*k*k2_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
uinc_re = d.Expression(cpp_code=cpp_g_re,
                       A=0.,                     # Will be overwritten later
                       k=0.,                     # Will be overwritten later
                       k0_=inc_dir[0],
                       k1_=inc_dir[1],
                       k2_=inc_dir[2],
                       degree=4)

cpp_uinc_im = ("A*k*k0_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
               "A*k*k1_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
               "A*k*k2_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
uinc_im = d.Expression(cpp_code=cpp_g_im,
                       A=0.,              # Will be overwritten later
                       k=0.,                   # Will be overwritten later
                       k0_=inc_dir[0],
                       k1_=inc_dir[1],
                       k2_=inc_dir[2],
                       degree=4)

# Spatial coordinate object to define piecewise functions within PML
x = d.SpatialCoordinate(mesh)

# Allocate loop-static pointer
#   for matrix A: discrete form of the bi-linear form and,
#   for vector b: discrete form of the linear form
A_mat = d.PETScMatrix()
b_vec = d.PETScVector()

# -------------------------------------------------------------------------------------------------
# LOOP OVER FREQUENCY
for frq_count, frq in enumerate(freqs):
    print('[*] Running solver for %d Hz' % frq)
    print('[*] Time : ', strftime("%Y-%m-%d %H:%M:%S", localtime()))
    h5counter = freqs_dict[frq]
    # ------------------------------------------------------
    # *** PARAMETERS for the given frequency *** #
    # Angular frequency
    omega = 2 * d.DOLFIN_PI * frq
    omega2 = omega ** 2
    # Wave number
    k = omega / c
    k2 = k ** 2

    # ------------------------------------------------------
    # *** VARIATIONAL FORM *** #

    # Define variational unknowns
    (u_re, u_im) = d.TrialFunctions(RTComplex)
    (v_re, v_im) = d.TestFunctions(RTComplex)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BOUNDARY CONDITIONS

    # # Apply Dirichlet boundary conditions
    # BC parameters
    # g = - u_inc :dot: n
    g_re.A = g_im.A = -1.0 / (rho * omega**2)
    g_re.k = g_im.k = k

    # # Dirichlet BCs to the functional space V
    # The only DBC is at object Boundary
    bcs = [d.DirichletBC(RTComplex.sub(0), g_re, boundaries, bd_marker["object"]),
           d.DirichletBC(RTComplex.sub(1), g_im, boundaries, bd_marker["object"]),
           d.DirichletBC(RTComplex.sub(0), zero_vec, boundaries, bd_marker["pml"]),
           d.DirichletBC(RTComplex.sub(1), zero_vec, boundaries, bd_marker["pml"])]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Define bi-linear form on LHS

    # # FLUID SUB-DOMAIN : Marker 0 # #
    # Within the bulk of Fluid domain
    # Bi-linear form
    a = (rho * c**2 * d.div(u_re) * d.div(v_re)
         + rho * c**2 * d.div(u_im) * d.div(v_im)
         - rho * omega2 * d.inner(u_re, v_re)
         - rho * omega2 * d.inner(u_im, v_im)) * dx(sd_marker["fluid"])

    # At the MPP Interface
    # MPP Impedance
    Z_mpp = mpp_impedance(omega)
    
    # Z_eff = i*Z_mpp
    Z_eff = C.Complex(float(Z_mpp.imag), -float(Z_mpp.real))
    
    # u :dot: n
    u_dot_n = C.Complex(d.inner(u_re('+'), n_vec('+')), d.inner(u_im('+'), n_vec('+')))
    
    # Bi-linear form at the MPP interface : Boundary Marker 2
    a += (omega * C.prod_re(Z_eff, u_dot_n) * d.inner(v_re('+'), n_vec('+')) * dS
          + omega * C.prod_im(Z_eff, u_dot_n) * d.inner(v_im('+'), n_vec('+')) * dS)

    # # PML SUB-DOMAIN : Marker 1 # #
    # PML absorption function
    sigma = lambda j: d.conditional(d.gt(abs(x[j]), L[j] + tol_s), c / abs(abs(x[j]) - (L[j] + Lpml)) / omega, zero)
    gamma = lambda j: C.Complex(d.Constant("1.0"), sigma(j))

    # Complex grad_U (only diagonal terms are relevant here)
    du_dx = lambda i: C.Complex(d.Dx(u_re[i], i), d.Dx(u_im[i], i))

    # Define Stiffness tensor * grad_U
    # S = C{i,j,k,l} * grad_U{k,l}
    S_re = sum(C.prod_re(C.Complex(1., 0.) / gamma(i), du_dx(i)) for i in range(3))
    S_im = sum(C.prod_im(C.Complex(1., 0.) / gamma(i), du_dx(i)) for i in range(3))

    # Define Mass Matrix * U
    # M_u = M * u
    u_j = lambda j: C.Complex(u_re[j], u_im[j])
    M_u_re = d.as_vector([C.prod_re(gamma(i), u_j(i)) for i in range(3)])
    M_u_im = d.as_vector([C.prod_im(gamma(i), u_j(i)) for i in range(3)])

    # Bi-linear form in the PML layer
    a += ((rho * c**2 * S_re * d.div(v_re)
           + rho * c**2 * S_im * d.div(v_im)
           - rho * omega2 * d.inner(M_u_re, v_re)
           - rho * omega2 * d.inner(M_u_im, v_im))
          * dx(sd_marker["pml"], scheme='default', degree=6))

    # # POROUS SUB-DOMAIN : Marker 2 # #
    # POROUS LAYER : Marker 2

    # Complex u, Div_U
    u = C.Complex(u_re, u_im)
    div_u = C.Complex(d.div(u_re), d.div(u_im))

    # Compute complex porous density and bulk modulus
    Rho = porous.rho(omega)
    BulkMod = porous.bulk_modulus(omega)
    CRho = C.Complex(float(Rho.real), float(Rho.imag))
    CBulkMod = C.Complex(float(BulkMod.real), float(BulkMod.imag))

    # Bi-linear form in the Porous layer
    a += (C.prod_re(CBulkMod, div_u) * d.div(v_re)
          + C.prod_im(CBulkMod, div_u) * d.div(v_im)
          - omega2 * d.inner(C.prod_re(CRho, u), v_re)
          - omega2 * d.inner(C.prod_im(CRho, u), v_im)) * dx(sd_marker["porous"])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Define linear form on RHS

    # Porous sub-domain
    factor_f_inc_re = float(Rho.real) * omega2 - float(BulkMod.real) * k2
    factor_f_inc_im = float(Rho.imag) * omega2 - float(BulkMod.imag) * k2
    factor_f_inc = C.Complex(factor_f_inc_re, factor_f_inc_im)

    uinc_re.A = uinc_im.A = inc_mag / (rho * omega2)
    uinc_re.k = uinc_im.k = k
    uinc = C.Complex(uinc_re, uinc_im)

    f_inc_re = d.as_vector(C.prod_re(factor_f_inc, uinc))
    f_inc_im = d.as_vector(C.prod_im(factor_f_inc, uinc))

    rhs = (d.inner(f_inc_re, v_re) + d.inner(f_inc_im, v_im)) * dx(sd_marker["porous"])

    # # At MPP-Interface
    u_inc_re_dot_n = d.inner(d.as_vector(g_re('+')), n_vec('+'))
    u_inc_im_dot_n = d.inner(d.as_vector(g_im('+')), n_vec('+'))
    u_inc_dot_n = C.Complex(u_inc_re_dot_n, u_inc_im_dot_n)
    
    rhs += (omega * C.prod_re(Z_eff, u_inc_dot_n) * d.inner(v_re('+'), n_vec('+')) * dS
           + omega * C.prod_im(Z_eff, u_inc_dot_n) * d.inner(v_im('+'), n_vec('+')) * dS)

    # rhs = (d.inner(zero_vec, v_re) * dx + d.inner(zero_vec, v_im) * dx)

    # ------------------------------------------------------
    # # SOLVE VARIATIONAL EQUATION

    # # Perform Matrix assembly on the variational form
    A = d.assemble(a, tensor=A_mat)
    b = d.assemble(rhs, tensor=b_vec)

    # Apply DirichletBCs to the assembled matrix
    for bc in bcs:
        bc.apply(A)
        bc.apply(b)

    # # Object to store solution
    w = d.Function(RTComplex)

    # # Invoke LASolver
    solver = d.LUSolver(A, 'mumps')
    solver.solve(w.vector(), b)

    # # Split real and imaginary parts of solution
    (u_sc_re, u_sc_im) = w.split(True)

    # # Name functions for reference and write to file
    u_sc_re.rename("Re(u_sc)", "Real Scattered field")
    u_sc_im.rename("Im(u_sc)", "Imag Scattered field")

    # ------------------------------------------------------
    # TOTAL FIELD COMPUTATIONS
    u_pw_re = uinc_re
    u_pw_im = uinc_im

    u_tot_re = d.interpolate(u_pw_re, RT)
    u_tot_im = d.interpolate(u_pw_im, RT)

    u_tot_re.vector().axpy(1.0, u_sc_re.vector())
    u_tot_im.vector().axpy(1.0, u_sc_im.vector())

    u_tot_re.rename("Re(u_tot)", "Real Total field")
    u_tot_im.rename("Im(u_tot)", "Imag Total field")

    # # ------------------------------------------------------
    # Write displacement fields in HDF5 for post-processing
    h5out.write(u_sc_re, "u_sc_re", h5counter)
    h5out.write(u_sc_im, "u_sc_im", h5counter)
    h5out.write(u_tot_re, "u_tot_re", h5counter)
    h5out.write(u_tot_im, "u_tot_im", h5counter)

h5out.close()
