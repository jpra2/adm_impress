from fully_implicit.input.mesh.Mesh_Manager import MeshManager
from fully_implicit.jacobian.Assembly import Ass
from impress import FineScaleMesh as mesh
from fully_implicit.jacobian.symbolic_jacobian import symbolic_J as symbolic_J
# M=mesh('./mesh/5x5x5.msh')
M=MeshManager('mesh/5x5x5.msh',dim=3)
# M=MeshManager('mesh/27x27x27.msh',dim=3)
# M=MeshManager('mesh/10x10x10.msh',dim=3)
# import pdb; pdb.set_trace()
J=symbolic_J().J
# import pdb; pdb.set_trace()
obj=Ass(M)
F_J=obj.F_Jacobian.copy()
# internos=obj.internal_faces.copy()
