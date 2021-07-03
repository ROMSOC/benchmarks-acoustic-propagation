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

# ------------------------------------------------------------------- #
#  Planewave Scattering by a rigid Sphere 
#  using the classical quadratic PML absorption function 
#
#  NOTE: 
#  It can also used for approximating the strength of (frequency-specific)
#  PML absorption coefficients, since exact solution is known and errors 
#  can be computed. Although, the monopole test case is suggested for it.
# ------------------------------------------------------------------- #
import dolfin as d
from complex import *
import os.path as path
import numpy as np

#d.set_log_level(d.LogLevel.DEBUG)
# -------------------------------------------------------------------------------------------------
# PARAMETERS
print('[*] Reading Parameters')

# PHYSICAL PARAMS
# Acoustic Fluid Properties
c = 340.  # sound speed
rho = 1.2

# MESH INPUT
meshfile = '../../00_meshes/sphere.xml.gz'

# MODEL PARAMS
# # Incident plane wave params
inc_mag = 1.  # Incident plane wave amplitude
inc_dir = np.array([1., 0., 0.])  # Incident plane wave direction

# # Geometry parameters
radius = 0.05
Lx = Ly = Lz = Lpml = 0.2

# PROBLEM PARAMETERS
freqs = np.array([800, 1200])  # np.arange(50, 401, 50)
nfreqs = freqs.size
omegas = 2 * np.pi * freqs

# PML calibration params
sigma0s = np.array([1.32e4, 1.05e5])

# IO PARAMS
# Directory where results are written. End with '/'
results_path = path.abspath('results/')

# -------------------------------------------------------------------------------------------------
# PROBLEM SETUP
print('[*] Setting up objects acc to Parameters')

# MESH
mesh = d.Mesh(meshfile)

# SUB-DOMAINS

# # Extract Sub-domain markers
subdomains = d.MeshFunction("size_t", mesh, mesh.topology().dim())

# # Sub-domain specifications
tol = 1e-10
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
pml_bnd = d.CompiledSubDomain("on_boundary && (fabs(x[0])>Lx+Lpml-tol || fabs(x[1])>Ly+Lpml-tol || fabs(x[2])>Lz+Lpml-tol)",
                              Lx=Lx, Ly=Ly, Lz=Lz, Lpml=Lpml, tol=tol)

object_bnd = d.CompiledSubDomain("on_boundary && fabs(x[0])<Lx+tol && fabs(x[1])<Ly+tol && fabs(x[2])<Lz+tol",
                                 Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)

# # Enumerate boundaries
boundary_markers.set_all(0)
object_bnd.mark(boundary_markers, 2)
pml_bnd.mark(boundary_markers, 1)

# # Obtain measures for boundary for boundary integrals
ds = d.Measure('ds', domain=mesh, subdomain_data=boundary_markers)

# NORMALS

# Compute normal vectors
n_vec = d.FacetNormal(mesh)

# FUNCTION SPACE SPECIFICATIONS

# # Define function space (Lagrange 1st polynomials)
P1 = d.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
Q = d.FunctionSpace(mesh, P1)

# # Define 2D-Real function space to substitute for Complex space
V = d.FunctionSpace(mesh, P1 * P1)   # NOTE: V.sub(0) and V.sub(1) are the subspaces

# INIT OUTPUT FILES
xdmf_file = d.XDMFFile(path.join(results_path, 'results.xdmf'))
xdmf_file.parameters['rewrite_function_mesh'] = False
xdmf_file.parameters['functions_share_mesh'] = True

# -------------------------------------------------------------------------------------------------
# DIAGNOSTIC OUTPUTS
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
print('[*] Setting up solver objects and expressions')

# Zero Expression (needed?)
zero = d.Constant("0.")

# Declare Variational form expressions
# # PML Absorption function
sx = d.Expression('fabs(x[0]) > Lx ? s0*pow(fabs(x[0])-Lx,2)/(Lpml*Lpml*w) : 0.',
                  Lx=Lx,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  degree=2)
sy = d.Expression('fabs(x[1]) > Ly ? s0*pow(fabs(x[1])-Ly,2)/(Lpml*Lpml*w) : 0.',
                  Ly=Ly,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  degree=2)
sz = d.Expression('fabs(x[2]) > Lz ? s0*pow(fabs(x[2])-Lz,2)/(Lpml*Lpml*w) : 0.',
                  Lz=Lz,
                  Lpml=Lpml,
                  w=0.,                      # Will be overwritten later
                  s0=0.,                     # Will be overwritten later
                  degree=2)

# Neumann boundary conditions
# # At Object Boundary
cpp_g_re = ("A*k*k0_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k1_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k2_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_re = d.Expression(cpp_code=cpp_g_re,
                    A=inc_mag,
                    k=0.,                     # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
                    degree=4)

cpp_g_im = ("-A*k*k0_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k1_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k2_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_im = d.Expression(cpp_code=cpp_g_im,
                    A=inc_mag,
                    k=0.,                   # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
                    degree=4)

# Field due to Plane Wave : $ p = A * exp( |k| \hat{k} \cdot \vec{r}) $
u_pw_re = d.Expression("A*cos(k*( k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                       A=inc_mag,
                       k=0.,                   # Will be overwritten later
                       k0_=inc_dir[0],
                       k1_=inc_dir[1],
                       k2_=inc_dir[2],
                       degree=4)
u_pw_im = d.Expression("A*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
                       A=inc_mag,
                       k=0.,                   # Will be overwritten later
                       k0_=inc_dir[0],
                       k1_=inc_dir[1],
                       k2_=inc_dir[2],
                       degree=4)

# Exact Solution expressions
cpp_uex_re = open("exact/uex_re.cpp").read()
cpp_uex_im = open("exact/uex_im.cpp").read()

# # Real part
uex_re_exp = d.CompiledExpression(d.compile_cpp_code(cpp_uex_re).ScatteringExact_Re(), degree=4)
uex_re_exp.p0 = inc_mag
uex_re_exp.a = radius

# # Imaginary Part
uex_im_exp = d.CompiledExpression(d.compile_cpp_code(cpp_uex_im).ScatteringExact_Im(), degree=4)
uex_im_exp.p0 = inc_mag
uex_im_exp.a = radius

# Total Field expressions
u_tot_exp = d.Expression('f1+f2',
                         f1=d.Function(Q),          # Will be overwritten later
                         f2=d.Function(Q),          # Will be overwritten later
                         degree=4)

# Error expressions
error_exp = d.Expression('f1-f2',
                         f1=d.Function(Q),               # Will be overwritten later
                         f2=d.Function(Q),               # Will be overwritten later
                         degree=4)
# Store L2 Errors
errors = np.zeros(nfreqs)

# LOOP OVER FREQUENCY
for ifreq in range(nfreqs):
    print('[*] Running solver for '+str(freqs[ifreq])+'Hz')

    # PARAMETERS for the given frequency
    # # Angular frequency
    omega = omegas[ifreq]

    # # Wave number
    k = omega / c
    k2 = k ** 2

    # # PML Absorpion coefficient
    sigma0 = sigma0s[ifreq]

    # EXACT SOLUTION PARAMETERS
    u_pw_re.k = u_pw_im.k = k
    uex_re_exp.k = uex_im_exp.k = k

    # VARIATIONAL FORM Parameters
    # # Set values of PML absorption function
    sx.w = sy.w = sz.w = omega
    sx.s0 = sy.s0 = sz.s0 = sigma0

    # # Diagonal stiffness tensors (real and imaginary part)
    # FIXME Multiplication of Complex Numbers yet to be defined
    gamma_x = Complex(d.Constant("1.0"), sx)
    gamma_y = Complex(d.Constant("1.0"), sy)
    gamma_z = Complex(d.Constant("1.0"), sz)

    diag1 = Complex(prod_re(gamma_y, gamma_z), prod_im(gamma_y, gamma_z))
    diag2 = Complex(prod_re(gamma_x, gamma_z), prod_im(gamma_x, gamma_z))
    diag3 = Complex(prod_re(gamma_x, gamma_y), prod_im(gamma_x, gamma_y))

    C_re = d.as_matrix(((div_re(diag1, gamma_x), zero, zero),
                        (zero, div_re(diag2, gamma_y), zero),
                        (zero, zero, div_re(diag3, gamma_z))))
    C_im = d.as_matrix(((div_im(diag1, gamma_x), zero, zero),
                        (zero, div_im(diag2, gamma_y), zero),
                        (zero, zero, div_im(diag3, gamma_z))))
    C = Complex(C_re, C_im)

    # # Mass coefficients
    # # FIXME M  = gamma_x * gamma_y * gamma_z
    M = Complex(prod_re(diag3, gamma_z), prod_im(diag3, gamma_z))

    # # Neumann BC parameters
    g_re.k = g_im.k = k

    # # Dirichlet BCs to the functional space V
    bcs = [d.DirichletBC(V.sub(0), zero, boundary_markers, 1),   # Exact Real at PML boundary
           d.DirichletBC(V.sub(1), zero, boundary_markers, 1)]   # Exact Imaginary at PML boundary

    # VARIATIONAL EQUATION Specifications
    # # Define variational unknowns
    (u_re, u_im) = d.TrialFunctions(V)
    (v_re, v_im) = d.TestFunctions(V)
    u = Complex(u_re, u_im)
    grad_u = Complex(d.grad(u_re), d.grad(u_im))

    # # Define variational problem
    a = (d.inner(prod_re(C, grad_u), d.grad(v_re)) * dx
         + d.inner(prod_im(C, grad_u), d.grad(v_im)) * dx
         - k2 * prod_re(M, u) * v_re * dx
         - k2 * prod_im(M, u) * v_im * dx)

    rhs = d.inner(g_re, n_vec) * v_re * ds(2) + d.inner(g_im, n_vec) * v_im * ds(2)

    # SOLVE VARIATIONAL FORM
    # # Perform Matrix assembly on the variational form
    A = d.assemble(a)
    b = d.assemble(rhs)

    # # Apply Dirichlet boundary conditions to the Assembly
    for bc in bcs:
        bc.apply(A)
        bc.apply(b)

    # # SOLVE VARIATIONAL EQUATION
    # # Object to store solution
    w = d.Function(V)

    # # Instantiate and invoke Solver
    solver = d.LUSolver(A, 'mumps')
    solver.solve(w.vector(), b)

    # # Split real and imaginary parts to obtain Scatter field
    (u_sc_re, u_sc_im) = w.split(True)
    u_sc_re.rename("Re(u_sc)", "FE Real scattered field")
    u_sc_im.rename("Im(u_sc)", "FE Imag scattered field")

    # Compute total field : Plane wave + scatter field
    u_tot_exp.f1, u_tot_exp.f2 = u_sc_re, u_pw_re
    u_re = d.interpolate(u_tot_exp, Q)
    u_tot_exp.f1, u_tot_exp.f2 = u_sc_im, u_pw_im
    u_im = d.interpolate(u_tot_exp, Q)
    u_re.rename("Re(u)", "FE Real computed")
    u_im.rename("Im(u)", "FE Imag. computed")

    # COMPUTE EXACT SOLUTION
    print('[ ] [*] Computing Exact solution')

    uex_re_interp = d.interpolate(uex_re_exp, Q)
    uex_im_interp = d.interpolate(uex_im_exp, Q)
    uex_re_interp.rename("Re(u_ex)", "Real exact")
    uex_im_interp.rename("Im(u_ex)", "Imag. exact")

    # COMPUTE ERRORS
    print('[ ] [*] Computing Errors')

    # # Compute Error field
    error_exp.f1, error_exp.f2 = u_re, uex_re_exp
    error_re = d.interpolate(error_exp, Q)
    error_exp.f1, error_exp.f2 = u_im, uex_im_exp
    error_im = d.interpolate(error_exp, Q)
    error_re.rename("Error_Re(u)", "Error Real ")
    error_im.rename("Error_Im(u)", "Error Imag")

    # # Compute Percentage relative error (%)
    errors[ifreq] = (100 * np.sqrt(d.assemble((error_re * error_re + error_im * error_im)
                                              * dx(0)))
                     / np.sqrt(d.assemble((uex_re_exp * uex_re_exp + uex_im_exp * uex_im_exp)
                                          * dx(0))))
    print("[ ] [*] L2-relative error (%): ", errors[ifreq])

    # OUTPUT RESULTS to file
    print('[ ] [*] Writing results to file')
    xdmf_file.write(u_re, ifreq)
    xdmf_file.write(u_im, ifreq)
    xdmf_file.write(uex_re_interp, ifreq)
    xdmf_file.write(uex_im_interp, ifreq)
    xdmf_file.write(u_sc_re, ifreq)
    xdmf_file.write(u_sc_im, ifreq)
    xdmf_file.write(error_re, ifreq)
    xdmf_file.write(error_im, ifreq)

xdmf_file.close()
np.savetxt(path.join(results_path, 'ERRORS.log'),
           np.array([freqs, sigma0s, errors]).T,
           fmt='%4.4f', delimiter='        ', newline='\n',
           header='Frequency     PML_Sigma_0       L2_Error')
