from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, ElementsData, DataManager, SparseDataManager, RockData, BiphasicData, SimulationData, WellsData, AccumulativeBiphasicData, CurrentBiphasicData
from packs.preprocess import TpfaPreprocess
from packs.properties import PhisicalProperties
from copy import deepcopy
import numpy as np
import os
import pdb
import time

dty = [('prod_o', float), ('prod_w', float), ('wor', float),
       ('delta_t', float), ('dvpi', float), ('loop', int), ('global_identifier', int)]


def preprocessar(M, data_impress, wells):

    phisical_properties = PhisicalProperties()

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
    elements.insert('volumes_to_faces', M.volumes.bridge_adjacencies(elements.volumes, 3, 2), 'array')

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

    biphasic_data = BiphasicData()
    biphasic_data['saturation'] = M.saturation[:].flatten()

    simulation_data = SimulationData()
    simulation_data['nkga_internal_faces'] = phisical_properties.get_nkga(rock_data['keq_faces'][elements.internal_faces], geom['u_direction_internal_faces'], geom['areas'][elements.internal_faces])


    current_data = CurrentBiphasicData()
    current_data['current'] = np.zeros(1, dtype=dty)

    accumulate = AccumulativeBiphasicData()
    accumulate.create()
    accumulate.insert_data(current_data['current'])

    DataManager.export_all_datas_to_npz()
    SparseDataManager.export_all_datas()

    return wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate, phisical_properties

def carregar():

    phisical_properties = PhisicalProperties()

    wells2 = WellsData(load=True)
    elements = ElementsData(load=True)
    geom = GeometricData(load=True)
    rock_data = RockData(load=True)
    biphasic_data = BiphasicData(load=True)
    simulation_data = SimulationData(load=True)
    current_data = CurrentBiphasicData(load=True)
    accumulate = AccumulativeBiphasicData()
    accumulate.create(global_identifier=current_data['current']['global_identifier'][0])

    return wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate, phisical_properties
