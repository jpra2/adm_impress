# from fully_implicit.input.mesh.Mesh_Manager import MeshManager
# from fully_implicit.get_jacobian.Assembly import Ass
from impress import FineScaleMesh as mesh
from fully_implicit.jacobian.symbolic_jacobian import symbolic_J as J
M=mesh('./mesh/5x5x5.msh')
# M=MeshManager('mesh/6x6.msh',dim=2)
symbolic_J=J().J
# import pdb; pdb.set_trace()
# obj=Ass(M)
# F_J=obj.F_Jacobian.copy()
# internos=obj.internal_faces.copy()
