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

import dolfin as d
from complex import *
import os.path as path
import numpy as np

## Optimization options for the form compiler
# d.parameters["form_compiler"]["optimize"] = True
# d.parameters["form_compiler"]["cpp_optimize"] = True
# -------------------------------------------------------------------------------------------------
# PARAMETERS
d.MPI.init()

if d.MPI.rank(d.MPI.comm_world) == 0:
    print('Number of ranks : ', d.MPI.size(d.MPI.comm_world))
    print('[*] Reading Parameters')

# PHYSICAL PARAMS
# Acoustic Fluid Properties
c = 340.  # sound speed
rho = 1.2

# MESH INPUT
meshfile = '../00_meshes/sphere.xml.gz'

# MODEL PARAMS
# # Incident plane wave params
inc_mag = 1.  # Incident plane wave amplitude
inc_dir = np.array([1., 0., 0.])  # Incident plane wave direction

# # # Sphere params
radius = 0.05
Lx = Ly = Lz = Lpml = 0.2


# PROBLEM PARAMETERS
freqs = np.array([400])
#freqs = np.arange(50, 401, 50)
nfreqs = freqs.size
omegas = 2 * np.pi * freqs

# PML calibration params
sigma0s = np.zeros(freqs.size)

# IO PARAMS
# Directory where results are written. End with '/'
results_path = path.abspath('./results/')

# -------------------------------------------------------------------------------------------------
# PROBLEM SETUP
if d.MPI.rank(d.MPI.comm_world) == 0:
    print('[*] Setting up objects acc to parameters')

# MESH
mesh = d.Mesh(meshfile)

# SUB-DOMAINS

# # Extract Sub-domain markers
subdomains = d.MeshFunction("size_t", mesh, mesh.topology().dim())

# # Sub-domain specifications
tol = d.DOLFIN_EPS_LARGE
pml_x = d.CompiledSubDomain("fabs(x[0])>Lx-tol && fabs(x[1])<Ly+tol && fabs(x[2])<Lz+tol",
                            Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_y = d.CompiledSubDomain("fabs(x[0])<Lx+tol && fabs(x[1])>Ly-tol && fabs(x[2])<Lz+tol",
                            Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_z = d.CompiledSubDomain("fabs(x[0])<Lx+tol && fabs(x[1])<Ly+tol && fabs(x[2])>Lz-tol",
                            Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_xy = d.CompiledSubDomain("fabs(x[0])>Lx-tol && fabs(x[1])>Ly-tol && fabs(x[2])<Lz+tol",
                             Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_yz = d.CompiledSubDomain("fabs(x[0])<Lx+tol && fabs(x[1])>Ly-tol && fabs(x[2])>Lz-tol",
                             Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_zx = d.CompiledSubDomain("fabs(x[0])>Lx-tol && fabs(x[1])<Ly+tol && fabs(x[2])>Lz-tol",
                             Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
pml_xyz = d.CompiledSubDomain("fabs(x[0])>Lx-tol && fabs(x[1])>Ly-tol && fabs(x[2])>Lz-tol",
                              Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)

# # Mark Sub-domains
subdomains.set_all(0)
pml_x.mark(subdomains, 1)
pml_y.mark(subdomains, 2)
pml_z.mark(subdomains, 3)
pml_xy.mark(subdomains, 4)
pml_yz.mark(subdomains, 5)
pml_zx.mark(subdomains, 6)
pml_xyz.mark(subdomains, 7)

# # Obtain volume integral measures for each sub-domain for volume integrals
dx = d.Measure('dx', domain=mesh, subdomain_data=subdomains)

# BOUNDARIES

# # Get all face markers
boundary_markers = d.MeshFunction("size_t", mesh, mesh.topology().dim()-1)

# # Specify expressions for boundaries
pml_bnd = d.CompiledSubDomain("on_boundary", Lx_=Lx+Lpml, Ly_=Ly+Lpml, Lz_=Lz+Lpml, tol=1e-6)
sphere_bnd = d.CompiledSubDomain("on_boundary && fabs(x[0])<Lx && fabs(x[1])<Ly && fabs(x[2])<Lz",
                                 Lx=Lx, Ly=Ly, Lz=Lz)

# # Enumerate boundaries
boundary_markers.set_all(0)
pml_bnd.mark(boundary_markers, 1)
sphere_bnd.mark(boundary_markers, 2)

# # Obtain measures for boundary for boundary integrals
ds = d.Measure('ds', domain=mesh, subdomain_data=boundary_markers)

# NORMALS

# Compute normal vectors
n_vec = d.FacetNormal(mesh)

# FUNCTION SPACE SPECIFICATIONS

# # Define function space (Raviart-Thomas elements)
RT = d.FiniteElement("RT", mesh.ufl_cell(), 1)
Q = d.FunctionSpace(mesh, RT)

# # Define 2D-Real function space to substitute for Complex space
V = d.FunctionSpace(mesh, RT * RT)   # NOTE: V.sub(0) and V.sub(1) are the subspaces


# INIT OUTPUT FILES

xdmf_file = d.XDMFFile(path.join(results_path, 'results.xdmf'))
xdmf_file.parameters['rewrite_function_mesh'] = False   # Mesh will not be written for every function
xdmf_file.parameters['functions_share_mesh'] = True     # All functions are over same mesh

# DIAGNOSTIC OUTPUTS
if d.MPI.rank(d.MPI.comm_world) == 0:
    # Write Problem Configuration log
    finfo = open(path.join(results_path, 'INFO.log'), 'w')
    finfo.write('MESH INFO:\n')
    finfo.write('Mesh Name = '+meshfile)
    finfo.write('\nNumber of vertices in mesh = ' + str(mesh.coordinates().shape[0])+'\n')
    finfo.write('h_min = ' + str(mesh.hmin()) + '\n')
    finfo.write('h_max = ' + str(mesh.hmax()) + '\n')
    finfo.write('\nPHYSICAL PARAMETERS : \n')
    finfo.write('Lx = '+str(Lx)+'\nLy = '+str(Ly)+'\nLz = '+str(Lz)+'\nLpml = '+str(Lpml)+'\n')
    finfo.write('\nFLUID PARAMETERS:\n')
    finfo.write('c = '+str(c)+'\nrho = '+str(rho)+'\n')
    finfo.write('\nPROBLEM CONFIG:\n')
    finfo.write('Freqs = '+str(freqs)+'\nSigma0 = '+str(sigma0s)+'\n')
    finfo.close()

# Write to file boundary and domain markers
fboundaries = d.File(path.join(results_path, 'BoundaryMarkers.pvd'))
fboundaries << boundary_markers
fdomains = d.File(path.join(results_path, 'DomainMarkers.pvd'))
fdomains << subdomains

# -------------------------------------------------------------------------------------------------
# SOLVER
if d.MPI.rank(d.MPI.comm_world) == 0:
    print('[*] Setting up solver objects and expressions')

# Zero Expression (needed?)
zero = d.Constant("0.")
zero_vec = d.Constant((0.,0., 0.))

# Declare Variational form expressions
# # PML Absorption function
sx = d.Expression('fabs(x[0]) > Lx+tol ? s0*pow(fabs(x[0])-Lx,2)/(Lpml*Lpml*w) : 0.',
                  Lx=Lx,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  tol=tol,
                  degree=0)
sy = d.Expression('fabs(x[1]) > Ly+tol ? s0*pow(fabs(x[1])-Ly,2)/(Lpml*Lpml*w) : 0.',
                  Ly=Ly,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  tol=tol,
                  degree=0)
sz = d.Expression('fabs(x[2]) > Lz+tol ? s0*pow(fabs(x[2])-Lz,2)/(Lpml*Lpml*w) : 0.',
                  Lz=Lz,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  tol=tol,
                  degree=0)

# # Neumann boundary conditions - Also Exact solution expression in this case
# # At Object Boundary
g_re = d.Expression(("A*k*k0_*sin(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                     "A*k*k1_*sin(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                     "A*k*k2_*sin(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))"),
                    A=0.,                     # Will be overwritten later
                    k=0.,                     # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
                    degree=4)
g_im = d.Expression(("-A*k*k0_*cos(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                     "-A*k*k1_*cos(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                     "-A*k*k2_*cos(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))"),
                    A=0.,                      # Will be overwritten later
                    k=0.,                      # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
                    degree=4)

# # Store L2 Errors
errors = np.zeros(nfreqs)

# Pressure field computation objects

# # Initialize Function Space for pressure field
Q2 = d.FunctionSpace(mesh, "DG", 0)

# Pressure fieldExact Solution expressions
# # Real part
pex_re = d.Expression("A*cos(-k*( k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                      A=inc_mag,
                      k=0.,                   # Will be overwritten later
                      k0_=inc_dir[0],
                      k1_=inc_dir[1],
                      k2_=inc_dir[2],
                      degree=0)
# # Imaginary Part
pex_im = d.Expression("A*sin(-k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                      A=inc_mag,
                      k=0.,                   # Will be overwritten later
                      k0_=inc_dir[0],
                      k1_=inc_dir[1],
                      k2_=inc_dir[2],
                      degree=0)
# Error expressions
error_p_exp = d.Expression('f1-f2',
                           f1=d.Function(Q2),               # Will be overwritten later
                           f2=d.Function(Q2),               # Will be overwritten later
                           degree=4)
# -------------------------------------------------------------------------------------------------
# LOOP OVER FREQUENCY
for ifreq in range(nfreqs):
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[*] Running solver for '+str(freqs[ifreq])+'Hz')
    # ------------------------------------------------------
    # PARAMETERS for the given frequency
    # # Angular frequency
    omega = omegas[ifreq]

    # # Wave number
    k = omega / c
    k2 = k ** 2

    # # PML Absorpion coefficient
    sigma0 = sigma0s[ifreq]

    # ------------------------------------------------------
    # VARIATIONAL FORM Parameters

    # # Define variational unknowns
    (u_re, u_im) = d.TrialFunctions(V)
    (v_re, v_im) = d.TestFunctions(V)
    u = Complex(u_re, u_im)

    # # Set values of PML absorption function
    sx.w = sy.w = sz.w = omega
    sx.s0 = sy.s0 = sz.s0 = sigma0

    gamma_x = Complex(d.Constant("1.0"), sx)
    gamma_y = Complex(d.Constant("1.0"), sy)
    gamma_z = Complex(d.Constant("1.0"), sz)
    gamma = [gamma_x, gamma_y, gamma_z]

    # # Mass matrix
    M = gamma_x * gamma_y * gamma_z

    # # Stiffness Tensor times grad U :  S = C{i,j,k,l} * grad_U{k,l}
    # # NOTE: Real and Imaginary parts of S are defined separately

    # # # Diagonal terms of grad_U :
    du_dx = lambda ure, uim, i: Complex(d.Dx(ure[i], i), d.Dx(uim[i], i))

    # # # Diagonal terms of S_Re
    S_ii_Re = lambda ure, uim, i: sum(prod_re(M / (gamma[i] * gamma[j]), du_dx(ure, uim, j)) for j in range(3))

    S_re = d.as_matrix([[S_ii_Re(u_re, u_im, j) if i == j else zero for i in range(3)] for j in range(3)])

    # # # Diagonal terms of S_Im
    S_ii_Im = lambda ure, uim, i: sum(prod_im(M / (gamma[i] * gamma[j]), du_dx(ure, uim, j)) for j in range(3))

    S_im = d.as_matrix([[S_ii_Im(u_re, u_im, j) if i == j else zero for i in range(3)] for j in range(3)])

    # # Define grad_v matrix
    # # NOTE: - Needs to be a UFL matrix for compatibility
    # #       - Defined only the relevant terms of grad_v (diagonal) to reduce assembling time
    T = lambda u: d.as_matrix([[d.Dx(u[i], i) if i == j else zero for i in range(3)] for j in range(3)])

    # # Neumann BC parameters
    g_re.A = g_im.A = inc_mag/(rho * omega**2)
    g_re.k = g_im.k = k

    # # Dirichlet BCs to the functional space V
    bcs = [d.DirichletBC(V.sub(0), g_re, boundary_markers, 1),
           d.DirichletBC(V.sub(1), g_im, boundary_markers, 1),
           d.DirichletBC(V.sub(0), g_re, boundary_markers, 2),
           d.DirichletBC(V.sub(1), g_im, boundary_markers, 2)]

    # ------------------------------------------------------
    # VARIATIONAL EQUATION Specifications

    # # Define bilinear form
    a = (d.inner(S_re, T(v_re)) * dx
         + d.inner(S_im, T(v_im)) * dx
         - k2 * d.inner(prod_re(M, u), v_re) * dx
         - k2 * d.inner(prod_im(M, u), v_im) * dx)

    # # Define linear form
    # # NOTE: RHS = 0, Defined this way for compatibility
    rhs = d.inner(zero_vec, v_re) * ds + d.inner(zero_vec, v_im) * ds

    # ------------------------------------------------------
    # SOLVE VARIATIONAL FORM

    # # Perform Matrix assembly on the variational form
    A = d.assemble(a)
    b = d.assemble(rhs)

    # # Apply Dirichlet boundary conditions
    for bc in bcs:
        bc.apply(A)
        bc.apply(b)

    # # SOLVE VARIATIONAL EQUATION
    # # Object to store solution
    w = d.Function(V)

    # # Instantiate and invoke Solver
    solver = d.LUSolver(A, 'mumps')
    solver.solve(w.vector(), b)

    # # Split real and imaginary parts
    (u_re, u_im) = w.split(True)
    u_re.rename("u_re", "Disp Real")
    u_im.rename("u_im", "Disp Imag")

    # ------------------------------------------------------
    # COMPUTE EXACT SOLUTION
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Computing exact solution')

    uex_re = d.interpolate(g_re, Q)
    uex_im = d.interpolate(g_im, Q)
    uex_re.rename("uex_re", "Disp_exact Real")
    uex_im.rename("uex_im", "Disp_exact exact")

    # ------------------------------------------------------
    # COMPUTE ERRORS
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Computing Errors')

    # # Compute Error field
    error_re = d.Function(Q)
    error_im = d.Function(Q)

    error_re.assign(u_re)
    error_re.vector().axpy(-1.0, uex_re.vector())
    error_re.rename("u_re_error", "Error Disp Real ")

    error_im.assign(u_im)
    error_im.vector().axpy(-1.0, uex_im.vector())
    error_im.rename("u_im_error", "Error Disp Imag")

    # # Compute Percentage relative error (%)
    error_rel = (100 * np.sqrt(d.assemble(d.dot(error_re, error_re) * dx(0)
                                          + d.dot(error_im, error_im) * dx(0)))
                 / np.sqrt(d.assemble(d.dot(uex_re, uex_re) * dx(0)
                                      + d.dot(uex_im, uex_im) * dx(0))))
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print("[ ] [*] L2-relative error (%): ", error_rel)

    # ------------------------------------------------------
    # OUTPUT RESULTS to file
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Writing results to file')

    xdmf_file.write(u_re, ifreq)
    xdmf_file.write(u_im, ifreq)
    xdmf_file.write(uex_re, ifreq)
    xdmf_file.write(uex_im, ifreq)
    xdmf_file.write(error_re, ifreq)
    xdmf_file.write(error_im, ifreq)

    # ------------------------------------------------------
    # PRESSURE FIELD COMPUTATIONS
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Computing Pressure field')

    # # Project Displacement solution on the DG0 Functionspace
    p_re = d.project(-rho * (c ** 2) * d.div(u_re), Q2)
    p_re.rename("p_re", "Pressure real")

    p_im = d.project(-rho * (c ** 2) * d.div(u_im), Q2)
    p_im.rename("p_im", "Pressure imag")

    # # Pressure field exact
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Computing Pressure field Exact')

    pex_re.k = pex_im.k = k
    pex_re_interp = d.interpolate(pex_re, Q2)
    pex_im_interp = d.interpolate(pex_im, Q2)
    pex_re_interp.rename("p_ex_re", "Pressure exact Real")
    pex_im_interp.rename("p_ex_im", "Pressure exact Imag")

    # # Compute errors in pressure field
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Computing pressure field errors')

    error_p_exp.f1, error_p_exp.f2 = p_re, pex_re
    error_p_re = d.interpolate(error_p_exp, Q2)
    error_p_exp.f1, error_p_exp.f2 = p_im, pex_im
    error_p_im = d.interpolate(error_p_exp, Q2)
    error_p_re.rename("Error_Re(p)", "Error Pr Real ")
    error_p_im.rename("Error_Im(p)", "Error Pr Imag")

    # # Compute Percentage relative error (%)
    error_p_rel = (100 * np.sqrt(d.assemble((error_p_re * error_p_re + error_p_im * error_p_im)
                                            * dx(0)))
                   / np.sqrt(d.assemble((pex_re * pex_re + pex_im * pex_im) * dx(0))))

    if d.MPI.rank(d.MPI.comm_world) == 0:
        print("[ ] [*] L2-relative error in pressure (%): ", error_p_rel)

    # # Write outputs to file
    if d.MPI.rank(d.MPI.comm_world) == 0:
        print('[ ] [*] Writing pressure fields to file')

    xdmf_file.write(p_re, ifreq)
    xdmf_file.write(p_im, ifreq)
    xdmf_file.write(pex_re_interp, ifreq)
    xdmf_file.write(pex_im_interp, ifreq)
    xdmf_file.write(error_p_re, ifreq)
    xdmf_file.write(error_p_im, ifreq)

xdmf_file.close()
