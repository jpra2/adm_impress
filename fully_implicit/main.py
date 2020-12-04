from input.mesh.Mesh_Manager import MeshManager
from Get_Jacobian.Assembly import Ass
import sympy as sym
from impress import FineScaleMesh as mesh


M=mesh('/mesh/5x5x5.msh')
# import pdb; pdb.set_trace()
M=MeshManager('mesh/6x6.msh',dim=2)
mesh=Ass(M)
F_J=mesh.F_Jacobian.copy()
internos=mesh.internal_faces.copy()
import pdb; pdb.set_trace()
# Jacobian=Assembly.Assembly(M)
# M1.mb.write_file("output/arq.vtk")
# M.imprima(M,"output/arq")s
#"Get_Jacobian/Assembly.py"
