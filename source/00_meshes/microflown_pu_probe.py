#!/usr/bin/env python

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

### Microflown PU Probe WITH PML LAYER Geometry + Mesh               ###

# NOTE 
# - Geometry includes Domain with and without Porous Layer. Choose the right Domain for meshing.
# - To run from console:  salome -t python path/to/sphere.py  (WARNING! This hasn't been tested!)
# - If running in GUI, click File > Load Script and select this file.
#   it may take a while for it to load since this script computes mesh before loading the study.
#   Comment out the compute lines for responsive behavior

###
### This file is generated automatically by SALOME v9.2.1 with dump python functionality
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


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(0.5, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 0.5, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 0.5)
Cylinder_1 = geompy.MakeCylinderRH(6.3, 15.95)
Vertex_1 = geompy.MakeVertex(0, -7, -2.65)
Cylinder_2 = geompy.MakeCylinder(Vertex_1, OY, 7.5, 14)
Bottom_part = geompy.MakeCommonList([Cylinder_1, Cylinder_2], True)
Vertex_2 = geompy.MakeVertex(0, -7, 19.05)
Cylinder_3 = geompy.MakeCylinder(Vertex_2, OY, 7.1, 14)
Top_part = geompy.MakeCommonList([Cylinder_1, Cylinder_3], True)
Vertex_3 = geompy.MakeVertex(0, -3.15, 3)
Vertex_4 = geompy.MakeVertex(0, 3.15, 3)
Cylinder_4 = geompy.MakeCylinder(Vertex_3, OZ, 2.5, 12)
Cylinder_5 = geompy.MakeCylinder(Vertex_4, OZ, 2.5, 12)
Pillar_1 = geompy.MakeCutList(Cylinder_4, [Bottom_part, Top_part], True)
Pillar_2 = geompy.MakeCutList(Cylinder_5, [Bottom_part, Top_part], True)
MPP_Covered_1 = geompy.MakeCutList(Cylinder_1, [Bottom_part, Top_part, Pillar_1, Pillar_2], True)
Vertex_9 = geompy.MakeVertex(0, 0, -73)
Probe_cylinder = geompy.MakeCylinder(Vertex_9, OZ, 6.3, 73)
Probe = geompy.MakePartitionNonSelfIntersectedShape([Bottom_part, Top_part, Pillar_2, Pillar_1, Probe_cylinder], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0, True)
Probe_centered = geompy.MakeTranslation(Probe, 0, 0, 28.52)
Mesh_Covered = geompy.MakeTranslation(MPP_Covered_1, 0, 0, 28.52)

# Fluid Area 
Vertex_10 = geompy.MakeVertex(20, 20, 60)
Vertex_11 = geompy.MakeVertex(-20, -20, -60)
Box_1 = geompy.MakeBoxTwoPnt(Vertex_10, Vertex_11)
Fluid = geompy.MakeCutList(Box_1, [Probe_centered, Mesh_Covered], True)

# # PML Area
Vertex_12 = geompy.MakeVertex(-20, -20, 75)
PML_Top = geompy.MakeBoxTwoPnt(Vertex_12, Vertex_10)
Vertex_13 = geompy.MakeVertex(20, 20, -75)
PML_Bot = geompy.MakeBoxTwoPnt(Vertex_11, Vertex_13)
Vertex_14 = geompy.MakeVertex(35, -20, -60)
PML_Rit = geompy.MakeBoxTwoPnt(Vertex_14, Vertex_10)
Vertex_15 = geompy.MakeVertex(-35, 20, 60)
PML_Lft = geompy.MakeBoxTwoPnt(Vertex_11, Vertex_15)
Vertex_16 = geompy.MakeVertex(-20, 35, -60)
PML_Frn = geompy.MakeBoxTwoPnt(Vertex_10, Vertex_16)
Vertex_17 = geompy.MakeVertex(20, -35, 60)
PML_Bck = geompy.MakeBoxTwoPnt(Vertex_17, Vertex_11)
Vertex_18 = geompy.MakeVertex(35, 35, -60)
PML_C1 = geompy.MakeBoxTwoPnt(Vertex_10, Vertex_18)
Vertex_19 = geompy.MakeVertex(-35, -35, 60)
PML_C2 = geompy.MakeBoxTwoPnt(Vertex_19, Vertex_11)
Vertex_20 = geompy.MakeVertex(35, -20, -60)
PML_C3 = geompy.MakeBoxTwoPnt(Vertex_20, Vertex_17)
Vertex_21 = geompy.MakeVertex(-35, 20, 60)
PML_C4 = geompy.MakeBoxTwoPnt(Vertex_21, Vertex_16)
Vertex_22 = geompy.MakeVertex(-20, 35, 75)
PML_C5 = geompy.MakeBoxTwoPnt(Vertex_22, Vertex_10)
Vertex_23 = geompy.MakeVertex(35, -20, 75)
PML_C6 = geompy.MakeBoxTwoPnt(Vertex_23, Vertex_10)
PML_C7 = geompy.MakeBoxTwoPnt(Vertex_15, Vertex_12)
PML_C8 = geompy.MakeBoxTwoPnt(Vertex_17, Vertex_12)
Vertex_24 = geompy.MakeVertex(20, -35, -75)
PML_C9 = geompy.MakeBoxTwoPnt(Vertex_24, Vertex_11)
PML_C10 = geompy.MakeBoxTwoPnt(Vertex_13, Vertex_14)
PML_C11 = geompy.MakeBoxTwoPnt(Vertex_13, Vertex_16)
Vertex_25 = geompy.MakeVertex(-35, 20, -75)
PML_C12 = geompy.MakeBoxTwoPnt(Vertex_25, Vertex_11)
Vertex_26 = geompy.MakeVertex(35, 35, 75)
PML_CC1 = geompy.MakeBoxTwoPnt(Vertex_10, Vertex_26)
Vertex_27 = geompy.MakeVertex(-35, -35, 60)
PML_CC2 = geompy.MakeBoxTwoPnt(Vertex_27, Vertex_12)
Vertex_28 = geompy.MakeVertex(35, -20, 75)
PML_CC3 = geompy.MakeBoxTwoPnt(Vertex_28, Vertex_17)
Vertex_29 = geompy.MakeVertex(-20, 35, 75)
PML_CC4 = geompy.MakeBoxTwoPnt(Vertex_29, Vertex_21)
Vertex_30 = geompy.MakeVertex(-35, -35, -75)
PML_CC5 = geompy.MakeBoxTwoPnt(Vertex_30, Vertex_11)
Vertex_31 = geompy.MakeVertex(35, 35, -60)
PML_CC6 = geompy.MakeBoxTwoPnt(Vertex_31, Vertex_13)
Vertex_32 = geompy.MakeVertex(20, -35, -75)
PML_CC7 = geompy.MakeBoxTwoPnt(Vertex_32, Vertex_14)
Vertex_33 = geompy.MakeVertex(-35, 20, -75)
PML_CC8 = geompy.MakeBoxTwoPnt(Vertex_33, Vertex_16)
Domain = geompy.MakePartitionNonSelfIntersectedShape([Mesh_Covered, Fluid, PML_Top, PML_Bot, PML_Rit, PML_Lft, PML_Frn, PML_Bck, PML_C1, PML_C2, PML_C3, PML_C4, PML_C5, PML_C6, PML_C7, PML_C8, PML_C9, PML_C10, PML_C11, PML_C12, PML_CC1, PML_CC2, PML_CC3, PML_CC4, PML_CC5, PML_CC6, PML_CC7, PML_CC8], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0, True)

# Domain with porous_layer
Vertex_34 = geompy.MakeVertex(0, 0, -25.53)
Cylinder_9 = geompy.MakeCylinder(Vertex_34, OZ, 19.78, 70)
Porous_Layer_3 = geompy.MakeCutList(Cylinder_9, [Probe_centered, Mesh_Covered], True)
Cylinder_10 = geompy.MakeCylinder(Vertex_34, OZ, 14.8, 70)
Cylinder_11 = geompy.MakeCylinder(Vertex_34, OZ, 9.65, 70)
Porous_Layer_2 = geompy.MakeCutList(Cylinder_10, [Probe_centered, Mesh_Covered], True)
Porous_layer_1 = geompy.MakeCutList(Cylinder_11, [Probe_centered, Mesh_Covered], True)
Fluid_with_mesh1_with_porous2 = geompy.MakeCutList(Box_1, [Probe_centered, Mesh_Covered, Porous_Layer_2], True)
Domain_porous2 = geompy.MakePartition([Mesh_Covered, PML_Top, PML_Bot, PML_Rit, PML_Lft, PML_Frn, PML_Bck, PML_C1, PML_C2, PML_C3, PML_C4, PML_C5, PML_C6, PML_C7, PML_C8, PML_C9, PML_C10, PML_C11, PML_C12, PML_CC1, PML_CC2, PML_CC3, PML_CC4, PML_CC5, PML_CC6, PML_CC7, PML_CC8, Porous_Layer_3, Fluid_with_mesh1_with_porous2], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
geompy.addToStudy( Bottom_part, 'Bottom_part' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Cylinder_3, 'Cylinder_3' )
geompy.addToStudy( Top_part, 'Top_part' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Cylinder_4, 'Cylinder_4' )
geompy.addToStudy( Cylinder_5, 'Cylinder_5' )
geompy.addToStudy( Pillar_1, 'Pillar_1' )
geompy.addToStudy( Pillar_2, 'Pillar_2' )
geompy.addToStudy( MPP_Covered_1, 'MPP_Covered_1' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Probe_cylinder, 'Probe_cylinder' )
geompy.addToStudy( Probe, 'Probe' )
geompy.addToStudy( Probe_centered, 'Probe_centered' )
geompy.addToStudy( Mesh_Covered, 'Mesh_Covered' )
geompy.addToStudy( Vertex_10, 'Vertex_10' )
geompy.addToStudy( Vertex_11, 'Vertex_11' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Fluid, 'Fluid' )
geompy.addToStudy( Vertex_12, 'Vertex_12' )
geompy.addToStudy( PML_Top, 'PML_Top' )
geompy.addToStudy( Vertex_13, 'Vertex_13' )
geompy.addToStudy( PML_Bot, 'PML_Bot' )
geompy.addToStudy( Vertex_14, 'Vertex_14' )
geompy.addToStudy( PML_Rit, 'PML_Rit' )
geompy.addToStudy( Vertex_15, 'Vertex_15' )
geompy.addToStudy( PML_Lft, 'PML_Lft' )
geompy.addToStudy( Vertex_16, 'Vertex_16' )
geompy.addToStudy( PML_Frn, 'PML_Frn' )
geompy.addToStudy( Vertex_17, 'Vertex_17' )
geompy.addToStudy( PML_Bck, 'PML_Bck' )
geompy.addToStudy( Vertex_18, 'Vertex_18' )
geompy.addToStudy( PML_C1, 'PML_C1' )
geompy.addToStudy( Vertex_19, 'Vertex_19' )
geompy.addToStudy( PML_C2, 'PML_C2' )
geompy.addToStudy( Vertex_20, 'Vertex_20' )
geompy.addToStudy( PML_C3, 'PML_C3' )
geompy.addToStudy( Vertex_21, 'Vertex_21' )
geompy.addToStudy( PML_C4, 'PML_C4' )
geompy.addToStudy( Vertex_22, 'Vertex_22' )
geompy.addToStudy( PML_C5, 'PML_C5' )
geompy.addToStudy( Vertex_23, 'Vertex_23' )
geompy.addToStudy( PML_C6, 'PML_C6' )
geompy.addToStudy( PML_C7, 'PML_C7' )
geompy.addToStudy( PML_C8, 'PML_C8' )
geompy.addToStudy( Vertex_24, 'Vertex_24' )
geompy.addToStudy( PML_C9, 'PML_C9' )
geompy.addToStudy( PML_C10, 'PML_C10' )
geompy.addToStudy( PML_C11, 'PML_C11' )
geompy.addToStudy( Vertex_25, 'Vertex_25' )
geompy.addToStudy( PML_C12, 'PML_C12' )
geompy.addToStudy( Vertex_26, 'Vertex_26' )
geompy.addToStudy( PML_CC1, 'PML_CC1' )
geompy.addToStudy( Vertex_27, 'Vertex_27' )
geompy.addToStudy( PML_CC2, 'PML_CC2' )
geompy.addToStudy( Vertex_28, 'Vertex_28' )
geompy.addToStudy( Vertex_29, 'Vertex_29' )
geompy.addToStudy( PML_CC3, 'PML_CC3' )
geompy.addToStudy( PML_CC4, 'PML_CC4' )
geompy.addToStudy( Vertex_30, 'Vertex_30' )
geompy.addToStudy( PML_CC5, 'PML_CC5' )
geompy.addToStudy( Vertex_31, 'Vertex_31' )
geompy.addToStudy( Vertex_32, 'Vertex_32' )
geompy.addToStudy( PML_CC6, 'PML_CC6' )
geompy.addToStudy( PML_CC7, 'PML_CC7' )
geompy.addToStudy( Vertex_33, 'Vertex_33' )
geompy.addToStudy( PML_CC8, 'PML_CC8' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudy( Vertex_34, 'Vertex_34' )
geompy.addToStudy( Cylinder_9, 'Cylinder_9' )
geompy.addToStudy( Porous_Layer_3, 'Porous_Layer_3' )
geompy.addToStudy( Cylinder_10, 'Cylinder_10' )
geompy.addToStudy( Cylinder_11, 'Cylinder_11' )
geompy.addToStudy( Porous_Layer_2, 'Porous_Layer_2' )
geompy.addToStudy( Porous_layer_1, 'Porous_layer_1' )
geompy.addToStudy( Fluid_with_mesh1_with_porous2, 'Fluid_with_mesh1_with_porous2' )
geompy.addToStudy( Domain_porous2, 'Domain_porous2' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
smeshObj_1 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
Mesh_1 = smesh.Mesh(Domain)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 1.5 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetChordalError( 0.1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetMinSize( 0.01 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.UnsetLocalSizeOnEntry("fluid")
NETGEN_3D_Parameters_1.UnsetLocalSizeOnEntry("fluid")
NETGEN_3D_Parameters_1.SetMinSize( 0.01 )
NETGEN_3D_Parameters_1.UnsetLocalSizeOnEntry("Domain")

isDone = Mesh_1.Compute()
Mesh_1.Reorient2DBy3D( [ Mesh_1 ], Mesh_1, 1 )

try:
  Mesh_1.ExportUNV( r'microflown_pu_mpp_mesh.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
