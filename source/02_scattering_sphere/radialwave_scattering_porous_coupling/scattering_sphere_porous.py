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
#   Radial wave scattering by a sphere with a spherical porous layer
#
# ------------------------------------------------------------------- #
import dolfin as d
import complex as C
import os.path as path
import numpy as np
import porous_JCAL as porous
from exact import wave_amplitudes as amplitudes

# Optimization options for the form compiler
d.parameters["form_compiler"]["optimize"] = True
d.parameters["form_compiler"]["cpp_optimize"] = True

# -------------------------------------------------------------------------------------------------
# PARAMETERS
print('[*] Reading Parameters')

# PHYSICAL PARAMS
# Acoustic Fluid Properties
c = 343.  # sound speed
rho = 1.213

# MESH INPUT
meshfile = '../../00_meshes/sphere_porous.xml.gz'

# MODEL PARAMS
# # Monopole params
monopole_mag = 1.  #

# Sphere params
radius = 0.05
Lx = Ly = Lz = Lpml = 0.2

# Porous layer params
r_min = 2*radius
r_max = 2.25*radius

# Boundary Condition parameters
# # Incident plane wave params
inc_mag = 1.  # Incident plane wave amplitude
inc_dir = np.array([1., 0., 0.])  # Incident plane wave direction

# PROBLEM PARAMETERS
freqs = np.array([100, 200, 400, 800, 1000])
nfreqs = freqs.size
omegas = 2 * np.pi * freqs

# IO PARAMS
# Directory where results are written. End with '/'
results_path = path.abspath('./results/')

# -------------------------------------------------------------------------------------------------
# PROBLEM SETUP
print('[*] Setting up objects acc to parameters')

# MESH
mesh = d.Mesh(meshfile)

# SUB-DOMAINS

# # Extract Sub-domain markers
subdomains = d.MeshFunction("size_t", mesh, mesh.topology().dim())

# # Sub-domain specifications
tol = d.DOLFIN_EPS_LARGE
pml_layer = d.CompiledSubDomain("fabs(x[0])>Lx-tol || fabs(x[1])>Ly-tol || fabs(x[2])>Lz-tol",
                                 Lx=Lx, Ly=Ly, Lz=Lz, tol=tol)
porous_layer = d.CompiledSubDomain("x.norm() > r_min-tol && x.norm() < r_max + tol",
                                   r_min=r_min, r_max=r_max, tol=tol)

# # Mark Sub-domains
subdomains.set_all(0)               # Fluid Domain = 0
pml_layer.mark(subdomains, 1)       # PML_Marker = 1
porous_layer.mark(subdomains, 2)    # Porous_Marker = 2

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
# # Define Function Space for pressure field
Q2 = d.FunctionSpace(mesh, "DG", 0)

P1 = d.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
Q3 = d.FunctionSpace(mesh, P1)

# # Define 2D-Real function space to substitute for Complex space
V = d.FunctionSpace(mesh, RT * RT)   # NOTE: V.sub(0) and V.sub(1) are the subspaces


# INIT OUTPUT FILES

xdmf_file = d.XDMFFile(path.join(results_path, 'results.xdmf'))
xdmf_file.parameters['rewrite_function_mesh'] = False   # Mesh will not be written for every function
xdmf_file.parameters['functions_share_mesh'] = True     # All functions are over same mesh

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
finfo.write('Freqs = '+str(freqs)+'\n')
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
zero_vec = d.Constant((0., 0., 0.))

# Declare Variational form expressions
x = d.SpatialCoordinate(mesh)

# # At Object Boundary
cpp_g_re = ("-A*k*k0_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k1_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "-A*k*k2_*sin(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_re = d.Expression(cpp_code=cpp_g_re,
                    A=0.,                     # Will be overwritten later
                    k=0.,                     # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
                    degree=4)

cpp_g_im = ("A*k*k0_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k1_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))",
            "A*k*k2_*cos(k*(k0_*x[0] + k1_*x[1] + k2_*x[2]))")
g_im = d.Expression(cpp_code=cpp_g_im,
                    A=0.,                   # Will be overwritten later
                    k=0.,                   # Will be overwritten later
                    k0_=inc_dir[0],
                    k1_=inc_dir[1],
                    k2_=inc_dir[2],
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

# Pressure field computation objects
bulkModulus_re_exp = d.Expression("""[&] () { auto r = x.norm();
                                          if ((r<r3)&&(r>r2)) {    
                                            return bulk_mod_re ; 
                                          } else {
                                            return rho_f * c_f * c_f; }  
                                        } ()""",
                                  r2=r_min,
                                  r3=r_max,
                                  rho_f=rho,
                                  c_f=c,
                                  bulk_mod_re=0., # To be set later
                                  degree=4)
bulkModulus_im_exp = d.Expression("""[&] () { auto r = x.norm();
                                          if ((r<r3)&&(r>r2)) {    
                                            return bulk_mod_im ; 
                                          } else {
                                            return 0.0; }  
                                        } ()""",
                                  r2=r_min,
                                  r3=r_max,
                                  bulk_mod_im=0., # To be set later
                                  degree=4)


# Allocate memory for matrix A
A = d.PETScMatrix()

# Exact Solution Expressions
cpp_uex_re = open("exact/uex_re.cpp").read()
cpp_uex_im = open("exact/uex_im.cpp").read()

# Displacement field Exact Solution expressions
uex_re = d.CompiledExpression(d.compile_cpp_code(cpp_uex_re).exp_uex_re(), degree=4)
uex_im = d.CompiledExpression(d.compile_cpp_code(cpp_uex_im).exp_uex_im(), degree=4)
# # Set Fixed Parameters for Exact Expressions
uex_re.r1 = uex_im.r1 = radius
uex_re.r2 = uex_im.r2 = r_min
uex_re.r3 = uex_im.r3 = r_max

# -------------------------------------------------------------------------------------------------
# LOOP OVER FREQUENCY
for ifreq in range(nfreqs):
#ifreq = 0
    # ------------------------------------------------------
    # PARAMETERS for the given frequency
    # # Angular frequency
    omega = omegas[ifreq]

    # # Wave number
    k = omega / c
    k2 = k ** 2
    omega2 = omega ** 2

    print('[*] Running solver for ' + str(freqs[ifreq]) + 'Hz ')
    # ------------------------------------------------------
    # VARIATIONAL FORM Parameters

    # # Set values of PML absorption function
    sx = d.conditional(d.gt(abs(x[0]), Lx + tol), c / abs(abs(x[0]) - (Lx + Lpml)) / omega, zero)
    sy = d.conditional(d.gt(abs(x[1]), Ly + tol), c / abs(abs(x[1]) - (Ly + Lpml)) / omega, zero)
    sz = d.conditional(d.gt(abs(x[2]), Lz + tol), c / abs(abs(x[2]) - (Lz + Lpml)) / omega, zero)

    gamma_x = C.Complex(d.Constant("1.0"), sx)
    gamma_y = C.Complex(d.Constant("1.0"), sy)
    gamma_z = C.Complex(d.Constant("1.0"), sz)
    gamma = [gamma_x, gamma_y, gamma_z]

    # # Mass matrix
    M = gamma_x * gamma_y * gamma_z

    # # BC parameters
    g_re.A = g_im.A = -inc_mag / (rho * omega2)
    g_re.k = g_im.k = k

    # # Dirichlet BCs to the functional space V
    bcs = [d.DirichletBC(V.sub(0), zero_vec, boundary_markers, 1),    # BC at PML Boundary
           d.DirichletBC(V.sub(1), zero_vec, boundary_markers, 1),    # BC at PML Boundary
           d.DirichletBC(V.sub(0), g_re, boundary_markers, 2),        # BC at Object Boundary
           d.DirichletBC(V.sub(1), g_im, boundary_markers, 2)]        # BC at Object Boundary

    # ------------------------------------------------------
    # VARIATIONAL EQUATION Specifications

    # # Define variational unknowns
    (u_re, u_im) = d.TrialFunctions(V)
    (v_re, v_im) = d.TestFunctions(V)
    u = C.Complex(u_re, u_im)
    div_u = C.Complex(d.div(u_re), d.div(u_im))

    # FLUID LAYER  : Marker 0
    # Variational forms within the Fluid Layer
    a = (rho * c**2 * d.div(u_re) * d.div(v_re) * dx(0)
         + rho * c**2 * d.div(u_im) * d.div(v_im) * dx(0)
         - rho * omega2 * d.inner(u_re, v_re) * dx(0)
         - rho * omega2 * d.inner(u_im, v_im) * dx(0))

    # PML LAYER : Marker 1
    # # Define Stiffness Tensor times grad U :  S = C{i,j,k,l} * grad_U{k,l}
    # # NOTE: Real and Imaginary parts are defined separately
    # # # Diagonal terms of grad_U :
    du_dx = lambda ure, uim, i: C.Complex(d.Dx(ure[i], i), d.Dx(uim[i], i))

    # Divergence operator in the PML domain
    Div_re = sum(C.prod_re(C.Complex(1., 0.) / gamma[i], du_dx(u_re, u_im, i)) for i in range(3))
    Div_im = sum(C.prod_im(C.Complex(1., 0.) / gamma[i], du_dx(u_re, u_im, i)) for i in range(3))

    # Scaled PML displacement vector to be used in the mass matrix
    coef = lambda ure, uim, i: C.Complex(ure[i], uim[i])
    u_scaled_re = d.as_vector([C.prod_re(gamma[i], coef(u_re, u_im, i)) for i in range(3)])
    u_scaled_im = d.as_vector([C.prod_im(gamma[i], coef(u_re, u_im, i)) for i in range(3)])

    # # Define bilinear form in the PML layer
    a += (rho * c**2 * Div_re * d.div(v_re) * dx(1, scheme='default', degree=6)
          + rho * c**2 * Div_im * d.div(v_im) * dx(1, scheme='default', degree=6)
          - rho * omega2 * d.inner(u_scaled_re, v_re) * dx(1, scheme='default', degree=6)
          - rho * omega2 * d.inner(u_scaled_im, v_im) * dx(1, scheme='default', degree=6))

    # POROUS LAYER : Marker 2
    Rho = porous.rho(omega)
    BulkMod = porous.bulk_modulus(omega)
    CRho = C.Complex(float(Rho.real), float(Rho.imag))
    CBulkMod = C.Complex(float(BulkMod.real), float(BulkMod.imag))

    a += (C.prod_re(CBulkMod, div_u) * d.div(v_re) * dx(2)
          + C.prod_im(CBulkMod, div_u) * d.div(v_im) * dx(2)
          - omega2 * d.inner(C.prod_re(CRho, u), v_re) * dx(2)
          - omega2 * d.inner(C.prod_im(CRho, u), v_im) * dx(2))
    # ------------------------------------------------------
    # # Define linear form
    #
    factor_f_inc_re = float(Rho.real) * omega2 - float(BulkMod.real) * k2
    factor_f_inc_im = float(Rho.imag) * omega2 - float(BulkMod.imag) * k2
    factor_f_inc = C.Complex(factor_f_inc_re, factor_f_inc_im)

    uinc_re.A = uinc_im.A = inc_mag / (rho * omega2)
    uinc_re.k = uinc_im.k = k
    uinc = C.Complex(uinc_re, uinc_im)

    f_inc_re = d.as_vector(C.prod_re(factor_f_inc, uinc))
    f_inc_im = d.as_vector(C.prod_im(factor_f_inc, uinc))

    rhs = d.inner(f_inc_re, v_re) * dx(2) + d.inner(f_inc_im, v_im) * dx(2)
    # ------------------------------------------------------
    # SOLVE VARIATIONAL FORM

    # # Perform Matrix assembly on the variational form
    A = d.assemble(a, tensor=A)
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
    (u_sc_re, u_sc_im) = w.split(True)

    u_sc_re.rename("Re(u_sc)", "Real Scattered field")
    u_sc_im.rename("Im(u_sc)", "Imag Scattered field")

    xdmf_file.write(u_sc_re, ifreq)
    xdmf_file.write(u_sc_im, ifreq)

    # ------------------------------------------------------
    # COMPUTE EXACT SOLUTION
    print('[ ] [*] Computing exact solution')

    # Compute coefficients for exact solution - needs to be contiguous for C++ compatibility
    coeffs = amplitudes.compute(inc_mag, omega, c, rho, Rho, BulkMod, radius, r_min, r_max)
    # # Contiguity check
    if not coeffs.flags["C_CONTIGUOUS"]:
        coeffs = np.ascontiguousarray(coeffs)

    # Set the dependent Exact Solution Expression parameters
    uex_re.coeffs = uex_im.coeffs = coeffs
    uex_re.k_f = uex_im.k_f = k
    uex_re.k_p = uex_im.k_p = np.sqrt(Rho * omega2 / BulkMod)
    uex_re.rho_f = uex_im.rho_f = rho
    uex_re.rho_p = uex_im.rho_p = Rho
    uex_re.omega2 = uex_im.omega2 = omega2

    uex_re_interp = d.interpolate(uex_re, Q)
    uex_im_interp = d.interpolate(uex_im, Q)
    uex_re_interp.rename("Re(u_ex)", "Real exact")
    uex_im_interp.rename("Im(u_ex)", "Imag. exact")

    # ------------------------------------------------------
    # COMPUTE ERRORS
    print('[ ] [*] Computing Errors')

    # # Compute Error field
    error_re = d.Function(Q)
    error_im = d.Function(Q)

    error_re.assign(u_sc_re)
    error_re.vector().axpy(-1.0, uex_re_interp.vector())
    error_re.rename("Error_Re(u)", "Error Real ")

    error_im.assign(u_sc_im)
    error_im.vector().axpy(-1.0, uex_im_interp.vector())
    error_im.rename("Error_Im(u)", "Error Imag")

    # # Compute Percentage relative error (%)
    error_rel = (100 * np.sqrt(d.assemble(d.dot(error_re, error_re) * dx(0)
                                          + d.dot(error_im, error_im) * dx(0)
                                          + d.dot(error_re, error_re) * dx(2)
                                          + d.dot(error_im, error_im) * dx(2)))
                / np.sqrt(d.assemble(d.dot(uex_re, uex_re) * dx(0)
                                      + d.dot(uex_im, uex_im) * dx(0)
                                      + d.dot(uex_re, uex_re) * dx(2)
                                      + d.dot(uex_im, uex_im) * dx(2))))
    print("[ ] [*] L2-relative error (%): ", error_rel)
    #errors_disp[ifreq] = error_rel

    # ------------------------------------------------------
    # OUTPUT RESULTS to file
    print('[ ] [*] Writing results to file')

    xdmf_file.write(uex_re_interp, ifreq)
    xdmf_file.write(uex_im_interp, ifreq)
    xdmf_file.write(error_re, ifreq)
    xdmf_file.write(error_im, ifreq)
    # ------------------------------------------------------
    # PRESSURE FIELD COMPUTATIONS
    print('[ ] [*] Computing Pressure field')

    # # Project Displacement solution on the DG0 Functionspace
    bulkModulus_re_exp.bulk_mod_re = BulkMod.real
    bulkModulus_im_exp.bulk_mod_im = BulkMod.imag
    bulkModulus_exp = C.Complex(bulkModulus_re_exp, bulkModulus_im_exp)
    div_u_sc = C.Complex(d.div(u_sc_re), d.div(u_sc_im))

    p_sc_re = d.project(-C.prod_re(bulkModulus_exp, div_u_sc), Q2)
    p_sc_re.rename("p_sc_re", "Pressure Scattered real")

    p_sc_im = d.project(-C.prod_im(bulkModulus_exp, div_u_sc), Q2)
    p_sc_im.rename("p_sc_im", "Pressure Scattered imag")

    # # Write outputs to file
    print('[ ] [*] Writing pressure fields to file')

    xdmf_file.write(p_sc_re, ifreq)
    xdmf_file.write(p_sc_im, ifreq)

    # ---------------------------------------------------------
    p_sc_re2 = d.project(p_sc_re, Q3)
    p_sc_re2.rename("p_sc_re_CG1", "Pressure Scattered real2")

    p_sc_im2 = d.project(p_sc_im, Q3)
    p_sc_im2.rename("p_sc_im_CG1", "Pressure Scattered imag2")

    # # Write outputs to file
    print('[ ] [*] Writing pressure fields 2 to file')

    xdmf_file.write(p_sc_re2, ifreq)
    xdmf_file.write(p_sc_im2, ifreq)

xdmf_file.close()
