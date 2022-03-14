# -*- coding: utf-8 -*-

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


###    SPHERE w Spherical Porous layer and a cartesian PML LAYER  ###
###                        Geometry + Mesh                        ###

# NOTE
# 1) To run from console:  salome -t python path/to/sphere.py
# 2) If running in GUI, it may take a while for it to load since this script computes mesh before loading the study.
#     Comment out the compute lines for responsive behavior

### This file was generated automatically by SALOME v8.4.0 with dump python functionality and then modified to include parameters.
###

import salome

salome.salome_init()

import salome_notebook
notebook = salome_notebook.NoteBook()


###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

### PARAMETERS ###
# Sphere params
radius = 0.05

# Domain params
lx = ly = lz = lpml = 0.2

# Porous layer params
r_min = 2 * radius
r_max = 2.25 * radius

mesh_param_sph = 5      # In fractional units of Sphere Radius (Hmin), default=10
mesh_param_por = 3       # In fractional units of Porous layer thickness (rmax-rmin)
mesh_param_max = 2       # In fractional units of Lpml (Hmax)

### GEOMETRY ###
geompy = geomBuilder.New()
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
b0 = geompy.MakeVertex(lx, ly, lz)
b1 = geompy.MakeVertex(-lx, -ly, -lz)

# Spherical geometry
Sphere = geompy.MakeSphereR(radius)

# Porous Layer
porous_in = geompy.MakeSphereR(r_min)
porous_out = geompy.MakeSphereR(r_max)
porous_layer = geompy.MakeCutList(porous_out, [porous_in], True)

# Fluid Domain
B = geompy.MakeBoxTwoPnt(b0, b1)
Fluid1 = geompy.MakeCutList(B, [porous_out], True)
Fluid2 = geompy.MakeCutList(porous_in, [Sphere], True)

## PML layers
# X-direction subdomain (2)
be0 = b0
be1 = geompy.MakeVertex( lx+lpml, -ly, -lz)
bw0 = geompy.MakeVertex(-lx-lpml, ly, lz)
bw1 = b1

Be = geompy.MakeBoxTwoPnt(be0,be1)
Bw = geompy.MakeBoxTwoPnt(bw0,bw1)

# Y-direction subdomain (2)
bn0 = b0
bn1 = geompy.MakeVertex(-lx, ly+lpml, -lz)
bs0 = geompy.MakeVertex(lx, -ly-lpml, lz)
bs1 = b1

Bn = geompy.MakeBoxTwoPnt(bn0,bn1)
Bs = geompy.MakeBoxTwoPnt(bs0,bs1)

# Z-direction subdomain (2)
front_b0 = geompy.MakeVertex(-lx, -ly, lz+lpml)
front_b1 = b0
back_b0 = b1
back_b1 = geompy.MakeVertex(lx, ly, -lz-lpml)

front_B = geompy.MakeBoxTwoPnt(front_b0,front_b1)
back_B = geompy.MakeBoxTwoPnt(back_b0,back_b1)

# XY-direction subdomain (clockwise) (4)
bne0 = b0
bne1 = geompy.MakeVertex(lx+lpml,ly+lpml,-lz)

bse0 = bs0
bse1 = be1

bsw0 = geompy.MakeVertex(-lx-lpml,-ly-lpml,lz)
bsw1 = b1

bnw0 = bw0
bnw1 = bn1

Bne = geompy.MakeBoxTwoPnt(bne0,bne1)
Bse = geompy.MakeBoxTwoPnt(bse0,bse1)
Bsw = geompy.MakeBoxTwoPnt(bsw0,bsw1)
Bnw = geompy.MakeBoxTwoPnt(bnw0,bnw1)

# YZ-direction subdomain (clockwise) (4)
back_bn0 = bn1
back_bn1 = back_b1

back_bs0 = b1
back_bs1 = geompy.MakeVertex(lx,-ly-lpml,-lz-lpml)

front_bs0 = front_b0
front_bs1 = bs0

front_bn0 = geompy.MakeVertex(-lx,ly+lpml,lz+lpml)
front_bn1 = b0

back_Bn = geompy.MakeBoxTwoPnt(back_bn0,back_bn1)
back_Bs = geompy.MakeBoxTwoPnt(back_bs0,back_bs1)
front_Bs = geompy.MakeBoxTwoPnt(front_bs0,front_bs1)
front_Bn = geompy.MakeBoxTwoPnt(front_bn0,front_bn1)

# ZX-direction subdomain (clockwise) (4)
front_be0 = geompy.MakeVertex(lx+lpml,-ly,lz+lpml)
front_be1 = b0

front_bw0 = front_b0
front_bw1 = bw0

back_bw0 = b1
back_bw1 = geompy.MakeVertex(-lx-lpml,ly,-lz-lpml)

back_be0 = be1
back_be1 = back_b1

front_Be = geompy.MakeBoxTwoPnt(front_be0,front_be1)
front_Bw = geompy.MakeBoxTwoPnt(front_bw0,front_bw1)
back_Bw = geompy.MakeBoxTwoPnt(back_bw0,back_bw1)
back_Be = geompy.MakeBoxTwoPnt(back_be0,back_be1)

# XYZ-direction subdomain (8)
front_bne0 = geompy.MakeVertex(lx+lpml,ly+lpml,lz+lpml)
front_bne1 = b0

front_bse0 = front_be0
front_bse1 = bs0

front_bsw0 = front_b0
front_bsw1 = bsw0

front_bnw0 = front_bn0
front_bnw1 = bw0

back_bne0 = bne1
back_bne1 = back_b1

back_bse0 = be1
back_bse1 = back_bs1

back_bsw0 = b1
back_bsw1 = geompy.MakeVertex(-lx-lpml,-ly-lpml,-lz-lpml)

back_bnw0 = bn1
back_bnw1 = back_bw1


front_Bne = geompy.MakeBoxTwoPnt(front_bne0,front_bne1)
front_Bse = geompy.MakeBoxTwoPnt(front_bse0,front_bse1)
front_Bsw = geompy.MakeBoxTwoPnt(front_bsw0,front_bsw1)
front_Bnw = geompy.MakeBoxTwoPnt(front_bnw0,front_bnw1)
back_Bne = geompy.MakeBoxTwoPnt(back_bne0,back_bne1)
back_Bse = geompy.MakeBoxTwoPnt(back_bse0,back_bse1)
back_Bsw = geompy.MakeBoxTwoPnt(back_bsw0,back_bsw1)
back_Bnw = geompy.MakeBoxTwoPnt(back_bnw0,back_bnw1)

# Build Domain
# Make Compound and Glue common faces
# Domain_init = geompy.MakeCompound([Be, Bw, Bn, Bs, front_B, back_B, Bne, Bse,Bsw, Bnw, back_Bn, back_Bs, front_Bs, front_Bn, front_Be, front_Bw, back_Bw, back_Be, front_Bne, front_Bse, front_Bsw, front_Bnw, back_Bne, back_Bse, back_Bsw, back_Bnw, Fluid])
# Domain = geompy.MakeGlueFaces(Domain_init, 1e-05)
Domainlist = [Be, Bw, Bn, Bs, front_B, back_B, Bne, Bse,Bsw, Bnw,
              back_Bn, back_Bs, front_Bs, front_Bn, front_Be,
              front_Bw, back_Bw, back_Be, front_Bne, front_Bse,
              front_Bsw, front_Bnw, back_Bne, back_Bse, back_Bsw,
              back_Bnw, Fluid2, porous_layer, Fluid1]

Domain = geompy.MakePartition(Domainlist)

# Extract Relevant subshapes for making conformal meshes
FluidLayer1 = geompy.SubShapeSortedCentres(Domain, geompy.ShapeType["SOLID"], [15])
PorousLayer = geompy.SubShapeSortedCentres(Domain, geompy.ShapeType["SOLID"], [14])
SphereFace = geompy.SubShapeSortedCentres(FluidLayer1, geompy.ShapeType["FACE"], [1])
#print(SphereID)

# Add to Study (for viewing)
#Verts
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( b0, 'b0' )
geompy.addToStudy( b1, 'b1' )
geompy.addToStudy( be1, 'be1' )
geompy.addToStudy( bw0, 'bw0' )
geompy.addToStudy( bn1, 'bn1' )
geompy.addToStudy( bs0, 'bs0' )
geompy.addToStudy( front_b0, 'front_b0' )
geompy.addToStudy( back_b1, 'back_b1' )
geompy.addToStudy( bne1, 'bne1' )
geompy.addToStudy( bsw0, 'bsw0' )
geompy.addToStudy( front_be0, 'front_be0' )
geompy.addToStudy( back_bs1, 'back_bs1' )
geompy.addToStudy( front_bn0, 'front_bn0' )
geompy.addToStudy( back_bw1, 'back_bw1' )
geompy.addToStudy( front_bne0, 'front_bne0' )
geompy.addToStudy( back_bsw1, 'back_bsw1' )

# Domain Components
geompy.addToStudy( Sphere, 'Sphere' )
geompy.addToStudy( B, 'b' )
geompy.addToStudy( Be, 'be' )
geompy.addToStudy( Bw, 'bw' )
geompy.addToStudy( Bn, 'bn' )
geompy.addToStudy( Bs, 'bs' )
geompy.addToStudy( front_B, 'front_b' )
geompy.addToStudy( back_B, 'back_b' )
geompy.addToStudy( Bne, 'bne' )
geompy.addToStudy( Bse, 'bse' )
geompy.addToStudy( Bsw, 'bsw' )
geompy.addToStudy( Bnw, 'bnw' )
geompy.addToStudy( back_Bn, 'back_bn' )
geompy.addToStudy( back_Bs, 'back_bs' )
geompy.addToStudy( front_Bs, 'front_bs' )
geompy.addToStudy( front_Bn, 'front_bn' )
geompy.addToStudy( front_Be, 'front_be' )
geompy.addToStudy( front_Bw, 'front_bw' )
geompy.addToStudy( back_Bw, 'back_bw' )
geompy.addToStudy( back_Be, 'back_be' )
geompy.addToStudy( front_Bne, 'front_bne' )
geompy.addToStudy( front_Bse, 'front_bse' )
geompy.addToStudy( front_Bsw, 'front_bsw' )
geompy.addToStudy( front_Bnw, 'front_bnw' )
geompy.addToStudy( back_Bne, 'back_bne' )
geompy.addToStudy( back_Bse, 'back_bse' )
geompy.addToStudy( back_Bsw, 'back_bsw' )
geompy.addToStudy( back_Bnw, 'back_bnw' )
geompy.addToStudy( Fluid1, 'Fluid Domain1' )
geompy.addToStudy( Fluid2, 'Fluid Domain2' )
geompy.addToStudy( porous_layer, 'Porous Domain' )
#geompy.addToStudy( Domain_init, 'Domain init' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, FluidLayer1, 'FluidLayer1' )
geompy.addToStudyInFather( Domain, PorousLayer, 'PorousLayer' )
geompy.addToStudyInFather( FluidLayer1, SphereFace, 'SphereFace' )

###
### SMESH component
### Mesh the created geometry

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
DomainMesh = smesh.Mesh(Domain)
# Tetrahedral mesh with NETGEN algorithm
NETGEN_1D_2D_3D = DomainMesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()

NETGEN_3D_Parameters_1.SetMaxSize(lpml/mesh_param_max)
NETGEN_3D_Parameters_1.SetMinSize(radius/mesh_param_sph)
NETGEN_3D_Parameters_1.SetLocalSizeOnShape(SphereFace, radius/mesh_param_sph)
NETGEN_3D_Parameters_1.SetLocalSizeOnShape(PorousLayer, (r_max-r_min)/mesh_param_por)


NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 3 )
NETGEN_3D_Parameters_1.SetChordalError( 0.1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )


isDone = DomainMesh.Compute()

# Reorient Faces to point their normals outside the volume
DomainMesh.Reorient2DBy3D( [DomainMesh], DomainMesh, 1 )

# Print information about the mesh
print("Information about mesh:")
print("Number of nodes       : ", DomainMesh.NbNodes())
print("          edges       : ", DomainMesh.NbEdges())
print("          faces       : ", DomainMesh.NbFaces())
print("          tetrahedrons: ", DomainMesh.NbTetras())
#
# Export Mesh as UNV file
try:
    #DomainMesh.ExportUNV(path
    #                     + '/results/sphere_porous'
    #                     + str(mesh_param_sph)
    #                     + str(mesh_param_por)
    #                     + str(mesh_param_max)
    #                     + '.unv' )
    DomainMesh.ExportUNV('sphere_porous.unv')
    #DomainMesh.ExportMED(path + '/results/sphere_refined' + str(mesh_param_max) + '.med')
    pass
except:
    #print('ExportUNV() failed. Invalid path specified :', path + '/unv/sphere_porous' + str(mesh_param_max) +'.unv')
    print('ExportUNV() failed. Invalid path specified :', 'sphere_porous.unv')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D params with refine')
smesh.SetName(DomainMesh.GetMesh(), 'DomainMesh')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
