from input.mesh.Mesh_Manager import MeshManager
from jacobian.Assembly import Ass
import sympy as sym
from impress import FineScaleMesh as mesh

import pdb; pdb.set_trace()
M=mesh('/5x5x5.h5m')
# import pdb; pdb.set_trace()
# M=MeshManager('mesh/6x6.msh',dim=2)
mesh=Ass(M)
F_J=mesh.F_Jacobian.copy()
internos=mesh.internal_faces.copy()
import pdb; pdb.set_trace()
# Jacobian=Assembly.Assembly(M)
# M1.mb.write_file("output/arq.vtk")
# M.imprima(M,"output/arq")s
#"Get_Jacobian/Assembly.py"
