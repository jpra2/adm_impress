
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh

from .. import directories as direc
mesh_name = direc.data_loaded['mesh_name']


M = msh(mesh_name, dim = 3)
