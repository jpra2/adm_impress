from packs.data_class import GeometricData
from packs.load.preprocessor0 import M
import pdb
import numpy as np

geom = GeometricData()
geom['unit_normal_vector'] = M.faces.normal
geom['areas'] = np.repeat(1.0, len(M.faces.all))
geom['block_dimension'] = np.repeat(np.array([1.0, 1.0, 1.0]), len(M.volumes.all))
