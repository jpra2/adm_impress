from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, ElementsData, DataManager, SparseDataManager, RockData, BiphasicData, SimulationData, WellsData
from packs.preprocess import TpfaPreprocess
from packs.running.initial_mesh_properties import initial_mesh
from packs.properties import PhisicalProperties
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from packs.tpfa.monophasic import MonophasicFluidProperties
from copy import deepcopy
import numpy as np
import os
import pdb
import time

M, elements_lv0, data_impress, wells = initial_mesh()
meshset_volumes = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
meshset_faces = M.core.mb.create_meshset()

def print_test_volumes(file_name):
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file(file_name, [meshset_volumes])

###########################################
## preprocessamento
def preprocessar():

    wells2 = WellsData()
    wells2._data = deepcopy(wells._data)

    perms = M.permeability[:]
    perms = perms.reshape(len(perms), 3, 3)

    elements = ElementsData()
    elements.set_mesh_elements(
        M.volumes.all,
        M.faces.all,
        M.edges.all,
        M.nodes.all,
        M.faces.boundary,
        M.nodes.boundary
    )
    elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
    elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
    elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))
    elements.insert('volumes_adj_internal_faces', elements.faces_to_volumes(elements.internal_faces), 'array')
    elements.insert('volumes_adj_boundary_faces', elements.faces_to_volumes(elements.boundary_faces).flatten(), 'array')
    elements.insert('volumes_adj_volumes_by_faces', M.volumes.bridge_adjacencies(elements.volumes, 2, 3), 'array')
    map_internal_faces = np.repeat(-1, len(elements.faces))
    map_internal_faces[elements.internal_faces] = np.arange(len(elements.internal_faces))
    elements.insert('map_internal_faces', map_internal_faces, 'array')

    n_vols = len(elements.volumes)

    preprocess = TpfaPreprocess()
    geom = GeometricData()
    geom['abs_u_normal_faces'] = M.faces.normal[:]
    geom['block_dimension'] = np.repeat(1.0, n_vols*3).reshape(n_vols, 3)
    geom['volume'] = np.prod(geom['block_dimension'], axis=1)
    geom['centroid_volumes'] = M.volumes.center[:]
    geom['centroid_faces'] = M.faces.center[:]
    geom['centroid_edges'] = M.edges.center[:]
    geom['centroid_nodes'] = M.nodes.center[:]
    nodes_of_faces = elements.faces_to_nodes(elements.faces)
    points_faces = geom['centroid_nodes'][nodes_of_faces]
    geom['areas'] = preprocess.get_areas_faces(elements.faces, points_faces)
    geom['u_direction_internal_faces'] = preprocess.get_u_normal_internal_faces(elements.get('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces], geom['centroid_volumes'])
    geom['hi'] = preprocess.get_h_internal_faces(geom['block_dimension'], elements.get('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces])

    rock_data = RockData()
    rock_data['permeability'] = perms.copy()
    rock_data['permeability_volumes_internal_faces_direction'] = preprocess.get_k_volumes_internal_faces_direction(elements.get('volumes_adj_internal_faces'), perms, geom['abs_u_normal_faces'][elements.internal_faces])
    rock_data['porosity'] = M.poro[:].flatten()
    rock_data['keq_faces'] = preprocess.get_equivalent_permeability_faces_from_diagonal_permeability(
        elements.volumes,
        elements.faces,
        elements.get('volumes_adj_internal_faces'),
        elements.get('volumes_adj_boundary_faces'),
        geom['abs_u_normal_faces'],
        rock_data['permeability'],
        elements.internal_faces,
        geom['block_dimension']
    )

    DataManager.export_all_datas_to_npz()
    SparseDataManager.export_all_datas()

    return wells2, elements, geom, rock_data
#########################################################

########################################
## carregamento
def carregar():

    wells2 = WellsData(load=True)
    elements = ElementsData(load=True)
    geom = GeometricData(load=True)
    rock_data = RockData(load=True)

    return wells2, elements, geom, rock_data
########################################

wells2, elements, geom, rock_data = preprocessar()
# wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate = carregar()
# pdb.set_trace()
solver = SolverSp()

monophasic_fluid_properties = MonophasicFluidProperties()
tpfa_monophasic = TpfaMonophasic()

# pdb.set_trace()

v0 = M.volumes.bridge_adjacencies(M.volumes.all, 2, 3)
v1 = elements.get('volumes_adj_volumes_by_faces')

v00 = M.faces.bridge_adjacencies(M.faces.internal[:], 2, 3)
v11 = elements.get('volumes_adj_internal_faces')
v22 = elements.faces_to_volumes(elements.internal_faces)

# pdb.set_trace()

transmissibility_faces = tpfa_monophasic.get_transmissibility_faces(
    geom['areas'],
    elements.internal_faces,
    elements.boundary_faces,
    elements.get('volumes_adj_internal_faces'),
    elements.get('volumes_adj_boundary_faces'),
    monophasic_fluid_properties.mi,
    rock_data['keq_faces'],
    geom['hi']
)

t1 = transmissibility_faces

# pdb.set_trace()

transmissibility_faces2 = tpfa_monophasic.get_transmissibility_faces(
    geom['areas'],
    M.faces.internal[:],
    M.faces.boundary[:],
    M.faces.bridge_adjacencies(M.faces.internal[:], 2, 3),
    M.faces.bridge_adjacencies(M.faces.boundary[:], 2, 3).flatten(),
    monophasic_fluid_properties.mi,
    rock_data['keq_faces'],
    geom['hi']
)

t2 = transmissibility_faces2

test1 = np.allclose(t1, t2)
# pdb.set_trace()

T2 = tpfa_monophasic.mount_transmissibility_matrix(
    t1[M.faces.internal[:]],
    M.faces.internal[:],
    M.faces.bridge_adjacencies(M.faces.internal[:], 2, 3),
    M.volumes.all
)

T1 = tpfa_monophasic.mount_transmissibility_matrix(
    t1[elements.internal_faces],
    elements.internal_faces,
    elements.get('volumes_adj_internal_faces'),
    elements.volumes
)

g_source_total_volumes = np.zeros(len(elements.volumes))

T_with_boundary, b = tpfa_monophasic.get_linear_problem(
    wells2['ws_p'],
    wells2['ws_q'],
    wells2['values_p'],
    wells2['values_q'],
    g_source_total_volumes,
    T1
)

pressure = solver.direct_solver(T_with_boundary, b)

M.pressure[:] = pressure

M.core.print(folder='results', file='test1', extension=".vtk", config_input="input_cards/print_settings0.yml")




pdb.set_trace()
