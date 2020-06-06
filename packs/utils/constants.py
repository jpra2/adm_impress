import numpy as np
from ..directories import data_loaded

def init(M):
    global n_volumes
    global v0
    global internal_faces
    global n_internal_faces
    global g
    global z
    global pretransmissibility_internal_faces
    global Pf
    global porosity
    global Cf
    global Cw
    global Pw
    global Vbulk
    global R

    Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
    Cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
    Pw = np.array(data_loaded['compositional_data']['water_data']['Pw']).astype(float)
    Cw = np.array(data_loaded['compositional_data']['water_data']['Cw']).astype(float)
    R = 8.3144598
    n_volumes = len(M.volumes.all)
    v0 = M.faces.bridge_adjacencies(M.faces.internal,2,3)
    porosity = M.data['poro']
    Vbulk = M.data['volume']
    internal_faces = M.faces.internal
    n_internal_faces = len(v0[:,0])
    g = 9.80665
    z = -M.data['centroid_volumes'][:,2]
    pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
    pretransmissibility_internal_faces = pretransmissibility_faces[ M.faces.internal]#[100]*np.ones(len(self.internal_faces))
