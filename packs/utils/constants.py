import numpy as np
from ..directories import data_loaded
from ..compositional import equation_of_state

def init(M, wells):
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
    global Vbulk
    global R
    global EOS_class
    global MUSCL
    global FR
    global bhp_ind
    global vols_no_wells

    EOS_class = getattr(equation_of_state, data_loaded['compositional_data']['equation_of_state'])
    MUSCL = data_loaded['compositional_data']['MUSCL']['set']
    FR = data_loaded['compositional_data']['FR']['set']
    Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
    Cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
    R = 8.3144598
    n_volumes = len(M.volumes.all)
    v0 = M.faces.bridge_adjacencies(M.faces.internal,2,3)
    porosity = M.data['poro']
    Vbulk = M.data['volume']
    internal_faces = M.faces.internal
    n_internal_faces = len(v0[:,0])
    g = 9.80665
    z = -M.data['centroid_volumes'][:,2]
    vols_index = M.volumes.all
    vols_no_wells = np.setdiff1d(vols_index,wells['all_wells'])
    pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
    pretransmissibility_internal_faces = pretransmissibility_faces[M.faces.internal]#[100]*np.ones(len(self.internal_faces))
    if len(wells['ws_p'])>1:
        bhp_ind = np.argwhere(M.volumes.center[wells['ws_p']][:,2] ==
        min(M.volumes.center[wells['ws_p']][:,2])).ravel()
    else: bhp_ind = wells['ws_p']

def component_properties():
    global load_k
    global load_w
    global compressible_k
    global n_phases
    global w
    global Bin
    global Tc
    global Pc
    global vc
    global Mw
    global Nc
    global s
    global n_components
    global Mw_w
    global Cw
    global Pw

    load_k = data_loaded['hidrocarbon_components']
    load_w = data_loaded['water_component']
    compressible_k = data_loaded['compressible_fluid']
    n_phases = 2 * load_k + 1 * load_w

    if load_k:
        w = np.array(data_loaded['compositional_data']['component_data']['w']).astype(float)
        Bin = np.array(data_loaded['compositional_data']['component_data']['Bin']).astype(float)
        Tc = np.array(data_loaded['compositional_data']['component_data']['Tc']).astype(float)
        Pc = np.array(data_loaded['compositional_data']['component_data']['Pc']).astype(float)
        vc = np.array(data_loaded['compositional_data']['component_data']['vc']).astype(float)
        Mw = np.array(data_loaded['compositional_data']['component_data']['Mw']).astype(float)
        s = np.array(data_loaded['compositional_data']['component_data']['vshift_parameter']).astype(float)
        Nc = len(Mw)
        #Pb_guess = np.array(data_loaded['compositional_data']['component_data']['Pb_guess']).astype(float)

    else: Nc = 0; z = []
    if load_w:
        Pw = np.array(data_loaded['compositional_data']['water_data']['Pw']).astype(float)
        Cw = np.array(data_loaded['compositional_data']['water_data']['Cw']).astype(float)
        Mw_w = data_loaded['compositional_data']['water_data']['Mw_w'] #* np.ones(n_volumes)
    else: Cw = 0
    n_components = Nc + 1 * load_w
