import numpy as np
from ..directories import data_loaded
from ..compositional import equation_of_state, equation_of_state_teste

def init(M, wells):
    global n_volumes
    global v0
    global in_vols_pairs
    global Ns_int_faces
    global N_sign
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
    global RS
    global bhp_ind
    global vols_no_wells
    global ds_faces
    global P_SC
    global T_SC
    global n_points
    global time_integration
    global hyperbolic_method


    P_SC = 101325
    T_SC = 288.706
    if data_loaded['water_miscible']:
        EOS_class = getattr(equation_of_state_teste, data_loaded['compositional_data']['equation_of_state'])
    else: EOS_class = getattr(equation_of_state, data_loaded['compositional_data']['equation_of_state'])

    MUSCL = dict()
    MUSCL['set'] = data_loaded['compositional_data']['MUSCL']['set']
    if MUSCL['set']:
        if data_loaded['compositional_data']['MUSCL']['minmod']:
            MUSCL['lim'] = 'minmod'
        elif data_loaded['compositional_data']['MUSCL']['VanAlbada1']:
            MUSCL['lim'] = 'Van_Albada1'
        elif data_loaded['compositional_data']['MUSCL']['VanLeer']:
            MUSCL['lim'] = 'Van_Leer'
        else: raise NameError('Missing inform the limiting strategy for MUSCL')

    FR = data_loaded['compositional_data']['FR']['set']

    if data_loaded['compositional_data']['MUSCL']['set']:
        hyperbolic_method = 'MUSCL'
    elif data_loaded['compositional_data']['FR']['set']:
        hyperbolic_method = 'FR'

    if MUSCL['set']: n_points = 2
    elif FR: n_points = data_loaded['compositional_data']['FR']['order']
    else: n_points=1

    #RS = dict()
    if data_loaded['compositional_data']['RiemannSolver']['LLF']:
        RS = 'LLF'
    elif data_loaded['compositional_data']['RiemannSolver']['MDW']:
        RS = 'MDW'
    elif data_loaded['compositional_data']['RiemannSolver']['DW']:
        RS = 'DW'
    elif data_loaded['compositional_data']['RiemannSolver']['ROE']:
        RS = 'ROE'
    else:
        RS = 'UPW'
        print('No Riemann Solver was marked -> default: traditional upwind')

    if data_loaded['compositional_data']['time_integration']['Euler']:
        time_integration = 'Euler'

    elif data_loaded['compositional_data']['time_integration']['RK3']:
        time_integration = 'RK3'
    else: raise NameError('Missing inform the time integration method used')


    Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
    Cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
    R = 8.3144598
    n_volumes = len(M.volumes.all)

    v00 = M.faces.bridge_adjacencies(M.faces.internal,2,3)
    c_int = M.faces.center(M.faces.internal)
    c_vols = M.volumes.center(M.volumes.all)
    try:
        pos = (c_int[:,np.newaxis,:] - c_vols[v00]).sum(axis=2)
        v0 = np.copy(v00)
        v0[:,0] = v00[pos>=0]
        v0[:,1] = v00[pos<0]
    except:
        v0 = v00

    in_vols_pairs, Ns_int_faces, N_sign = correct_ctev0(M)

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
    ds_faces_axis = M.data['centroid_volumes'][v0[:,1],:] -  M.data['centroid_volumes'][v0[:,0],:]
    ds_faces = ds_faces_axis.sum(axis=-1)

    if len(wells['ws_p'])>1:
        bhp_ind = np.argwhere(M.volumes.center[wells['ws_p']][:,2] ==
            min(M.volumes.center[wells['ws_p']][:,2])).ravel()
    else: bhp_ind = wells['ws_p']


def component_properties():
    global load_k
    global load_w
    global compressible_k
    global miscible_w
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
    global Csi_W_def

    Csi_W_def = 55549.81959478213

    load_k = data_loaded['hidrocarbon_components']
    load_w = data_loaded['water_component']
    compressible_k = data_loaded['compressible_oil']
    miscible_w = data_loaded['water_miscible']
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
    n_components = Nc + 1 * load_w * (1-miscible_w)


def correct_ctev0(M):


    internal_faces = M.faces.internal[:]

    # Retrieve the points that form the components of the normal vectors.
    internal_faces_nodes = M.data['faces_nodes']

    I_idx = internal_faces_nodes[:, 0]
    J_idx = internal_faces_nodes[:, 1]
    K_idx = internal_faces_nodes[:, 2]

    I = M.nodes.coords[I_idx]
    J = M.nodes.coords[J_idx]
    K = M.nodes.coords[K_idx]
    #I[...,-1] = abs(I[...,-1])
    #J[...,-1] = abs(J[...,-1])
    #K[...,-1] = abs(K[...,-1])

    internal_faces = M.faces.internal[:]
    in_vols_pairs = M.faces.bridge_adjacencies(
        internal_faces, 2, 3)

    n_vols_pairs = len(internal_faces)
    internal_volumes_centers_flat = M.volumes.center[in_vols_pairs.flatten(
    )]
    internal_volumes_centers = internal_volumes_centers_flat.reshape((
        n_vols_pairs, 2, 3))

    '''LJ  = J-internal_volumes_centers[:, 0]

    # Set the normal vectors.
    Ns = 0.5 * np.cross(I-J,K-J)
    Ns_norm = np.linalg.norm(Ns, axis=1)
    Ns /= Ns_norm[:,np.newaxis]
    import pdb; pdb.set_trace()

    N_sign = np.sign(np.einsum("ij,ij->i", LJ, Ns))

    (in_vols_pairs[N_sign < 0, 0],
     in_vols_pairs[N_sign < 0, 1]) = (in_vols_pairs[N_sign < 0, 1],
                                           in_vols_pairs[N_sign < 0, 0])
    import pdb; pdb.set_trace()'''

    #teste meu - Duda - ver se pode funcionar

    Ns = M.data['faces_normals']
    #Ns_norm = np.linalg.norm(Ns, axis=1)
    #Ns /= Ns_norm[:,np.newaxis]

    LJ = M.faces.center(M.faces.internal) -internal_volumes_centers[:, 0]
    N_sign = np.sign(np.einsum("ij,ij->i", LJ, Ns[internal_faces]))
    (in_vols_pairs[N_sign < 0, 0],
     in_vols_pairs[N_sign < 0, 1]) = (in_vols_pairs[N_sign < 0, 1],
                                           in_vols_pairs[N_sign < 0, 0])

    return in_vols_pairs, Ns, N_sign
