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

# Optimization options for the form compiler
d.parameters["form_compiler"]["optimize"] = True
d.parameters["form_compiler"]["cpp_optimize"] = True

# -------------------------------------------------------------------------------------------------
# PARAMETERS
print('[*] Reading Parameters')

# PHYSICAL PARAMS
# Acoustic Fluid Properties
c = 340.  # sound speed
rho = 1.2

# MESH INPUT
meshfile = '../00_meshes/sphere.xml.gz'

# MODEL PARAMS
# # Monopole params
monopole_mag = 1.  #

# Sphere params
radius = 0.05
Lx = Ly = Lz = Lpml = 0.2

# PROBLEM PARAMETERS
freqs = np.array([800])
#freqs = np.arange(50, 401, 50)
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
zero_vec = d.Constant((0.,0., 0.))

x = d.SpatialCoordinate(mesh)


# # Neumann boundary conditions - Also Exact solution expression in this case
# # At Object Boundary
cpp_g_re = ("""[&] () { auto r = x.norm();
                        return A * x[0] * (-k*r*sin(k*r) - cos(k*r)) / pow(r,3);
                } () """,
            """[&] () { auto r = x.norm();
                        return A * x[1] * (-k*r*sin(k*r) - cos(k*r)) / pow(r,3);
             } ()""",
            """[&] () { auto r = x.norm();
                        return A * x[2] * (-k*r*sin(k*r) - cos(k*r)) / pow(r,3);
             } () """)
g_re = d.Expression(cpp_code=cpp_g_re,
                    A=0.,                     # Will be overwritten later
                    k=0.,                     # Will be overwritten later
                    degree=3)

cpp_g_im = ("""[&] () { auto r = x.norm();
                        return A * x[0] * (k*r*cos(k*r) - sin(k*r)) / pow(r,3);
                } () """,
            """[&] () { auto r = x.norm();
                        return A * x[1] * (k*r*cos(k*r) - sin(k*r)) / pow(r,3);
             } ()""",
            """[&] () { auto r = x.norm();
                        return A * x[2] * (k*r*cos(k*r) - sin(k*r)) / pow(r,3);
             } () """)

g_im = d.Expression(cpp_code=cpp_g_im,
                    A=0.,                   # Will be overwritten later
                    k=0.,                   # Will be overwritten later
                    degree=3)

# Pressure field computation objects

# # Initialize Function Space for pressure field
Q2 = d.FunctionSpace(mesh, "DG", 0)

# Pressure fieldExact Solution expressions
# # Real part
cpp_pex_re = """
            [&] () {
                auto r = x.norm();
                return A * cos(k*r) / r;
            } ()
        """
pex_re = d.Expression(cpp_code=cpp_pex_re,
                      A=monopole_mag,
                      k=0.,                   # Will be overwritten later
                      degree=0)
# # Imaginary Part
cpp_pex_im = """
            [&] () {
                auto r = x.norm();
                return A * sin(k*r) / r;
            } ()
        """
pex_im = d.Expression(cpp_code=cpp_pex_im,
                      A=monopole_mag,
                      k=0.,                   # Will be overwritten later
                      degree=0)
# Error expressions
error_p_exp = d.Expression('f1-f2',
                           f1=d.Function(Q2),               # Will be overwritten later
                           f2=d.Function(Q2),               # Will be overwritten later
                           degree=4)

# Store Errors
errors_disp = np.zeros(freqs.size)
errors_pr = np.zeros(freqs.size)
# -------------------------------------------------------------------------------------------------
# LOOP OVER FREQUENCY
for ifreq in range(nfreqs):
    # ------------------------------------------------------
    # PARAMETERS for the given frequency
    # # Angular frequency
    omega = omegas[ifreq]

    # # Wave number
    k = omega / c
    k2 = k ** 2

    print('[*] Running solver for ' + str(freqs[ifreq]) + 'Hz')
    # ------------------------------------------------------
    # VARIATIONAL FORM Parameters

    # # Set values of BC expressions
    g_re.A = g_im.A = monopole_mag / (rho * omega ** 2)
    g_re.k = g_im.k = k

    # # Dirichlet BCs to the functional space V
    bcs = [d.DirichletBC(V.sub(0), zero_vec, boundary_markers, 1),  # BC at PML Boundary
           d.DirichletBC(V.sub(1), zero_vec, boundary_markers, 1),  # BC at PML Boundary
           d.DirichletBC(V.sub(0), g_re, boundary_markers, 2),  # BC at Object Boundary
           d.DirichletBC(V.sub(1), g_im, boundary_markers, 2)]  # BC at Object Boundary

    # PML stretching functions
    sx = d.conditional(d.gt(abs(x[0]), Lx + tol), c / abs(abs(x[0]) - (Lx + Lpml)) / omega, zero)
    sy = d.conditional(d.gt(abs(x[1]), Ly + tol), c / abs(abs(x[1]) - (Ly + Lpml)) / omega, zero)
    sz = d.conditional(d.gt(abs(x[2]), Lz + tol), c / abs(abs(x[2]) - (Lz + Lpml)) / omega, zero)

    gamma_x = Complex(d.Constant("1.0"), sx)
    gamma_y = Complex(d.Constant("1.0"), sy)
    gamma_z = Complex(d.Constant("1.0"), sz)
    gamma = [gamma_x, gamma_y, gamma_z]

    # # Mass matrix
    M = gamma_x * gamma_y * gamma_z

    # ------------------------------------------------------
    # VARIATIONAL EQUATION Specifications

    # # Define variational unknowns
    (u_re, u_im) = d.TrialFunctions(V)
    (v_re, v_im) = d.TestFunctions(V)
    u = Complex(u_re, u_im)

    # # Define Stiffness Tensor times grad U :  S = C{i,j,k,l} * grad_U{k,l}
    # # NOTE: Real and Imaginary parts are defined separately

    # # # Diagonal terms of grad_U :
    du_dx = lambda ure, uim, i: Complex(d.Dx(ure[i], i), d.Dx(uim[i], i))


    # # # Diagonal terms of S_Re
    S_ii_Re = lambda ure, uim, i: sum(prod_re(M/(gamma[i]*gamma[j]), du_dx(ure, uim, j)) for j in range(3))

    S_re = d.as_matrix([[S_ii_Re(u_re, u_im, j) if i == j else zero for i in range(3)] for j in range(3)])

    # # # Diagonal terms of S_Im
    S_ii_Im = lambda ure, uim, i: sum(prod_im(M/(gamma[i]*gamma[j]), du_dx(ure, uim, j)) for j in range(3))

    S_im = d.as_matrix([[S_ii_Im(u_re, u_im, j) if i == j else zero for i in range(3)] for j in range(3)])

    # # Define grad_v matrix
    # # NOTE: - Needs to be a UFL matrix for compatibility
    # #       - Defined only the relevant terms of grad_v (diagonal) to reduce assembling time
    T = lambda u: d.as_matrix([[d.Dx(u[i], i) if i == j else zero for i in range(3)] for j in range(3)])

    # # Define bilinear form
    a = (d.inner(S_re, T(v_re)) * dx
         + d.inner(S_im, T(v_im)) * dx
         - k2 * d.inner(prod_re(M, u), v_re) * dx
         - k2 * d.inner(prod_im(M, u), v_im) * dx)

    # # Define linear form
    # # NOTE: RHS = 0, Defined this way for compatibility
    rhs = d.inner(zero_vec, v_re) * dx + d.inner(zero_vec, v_im) * dx

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
    u_re.rename("Re(u)", "FEM Real field")
    u_im.rename("Im(u)", "FEM Imag field")

    # ------------------------------------------------------
    # COMPUTE EXACT SOLUTION
    print('[ ] [*] Computing exact solution')

    #g_re.A = g_im.A = monopole_mag / (rho * omega ** 2)
    uex_re = d.interpolate(g_re, Q)
    uex_im = d.interpolate(g_im, Q)
    uex_re.rename("Re(u_ex)", "Real exact")
    uex_im.rename("Im(u_ex)", "Imag. exact")

    # ------------------------------------------------------
    # COMPUTE ERRORS
    print('[ ] [*] Computing Errors')

    # # Compute Error field
    error_re = d.Function(Q)
    error_im = d.Function(Q)

    error_re.assign(u_re)
    error_re.vector().axpy(-1.0, uex_re.vector())
    error_re.rename("Error_Re(u)", "Error Real ")

    error_im.assign(u_im)
    error_im.vector().axpy(-1.0, uex_im.vector())
    error_im.rename("Error_Im(u)", "Error Imag")

    # # Compute Percentage relative error (%)
    error_rel = (100 * np.sqrt(d.assemble(d.dot(error_re, error_re) * dx(0)
                                          + d.dot(error_im, error_im) * dx(0)))
                 / np.sqrt(d.assemble(d.dot(uex_re, uex_re) * dx(0)
                                      + d.dot(uex_im, uex_im) * dx(0))))
    print("[ ] [*] L2-relative error (%): ", error_rel)
    errors_disp[ifreq] = error_rel
    # ------------------------------------------------------
    # OUTPUT RESULTS to file
    print('[ ] [*] Writing results to file')

    xdmf_file.write(u_re, ifreq)
    xdmf_file.write(u_im, ifreq)
    xdmf_file.write(uex_re, ifreq)
    xdmf_file.write(uex_im, ifreq)
    xdmf_file.write(error_re, ifreq)
    xdmf_file.write(error_im, ifreq)

    # ------------------------------------------------------
    # PRESSURE FIELD COMPUTATIONS
    print('[ ] [*] Computing Pressure field')

    # # Project Displacement solution on the DG0 Functionspace
    p_re = d.project(-rho * (c ** 2) * d.div(u_re), Q2)
    p_re.rename("p_re", "Pressure real")

    p_im = d.project(-rho * (c ** 2) * d.div(u_im), Q2)
    p_im.rename("p_im", "Pressure imag")

    # # Pressure field exact
    print('[ ] [*] Computing Pressure field Exact')

    pex_re.k = pex_im.k = k
    pex_re_interp = d.interpolate(pex_re, Q2)
    pex_im_interp = d.interpolate(pex_im, Q2)
    pex_re_interp.rename("p_ex_re", "Pressure exact Real")
    pex_im_interp.rename("p_ex_im", "Pressure exact Imag")

    # # Compute errors in pressure field
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

    print("[ ] [*] L2-relative error in pressure (%): ", error_p_rel)
    errors_pr[ifreq] = error_p_rel

    # # Write outputs to file
    print('[ ] [*] Writing pressure fields to file')

    xdmf_file.write(p_re, ifreq)
    xdmf_file.write(p_im, ifreq)
    xdmf_file.write(pex_re_interp, ifreq)
    xdmf_file.write(pex_im_interp, ifreq)
    xdmf_file.write(error_p_re, ifreq)
    xdmf_file.write(error_p_im, ifreq)

xdmf_file.close()
np.savetxt(path.join(results_path, 'ERRORS.log'),
           np.array([freqs, errors_disp, errors_pr]).T,
           fmt='%4.4f', delimiter='        ', newline='\n',
           header='Frequency     L2_Error_Disp  L2_Error_pr')