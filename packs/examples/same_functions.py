from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.generic_data import PrimalCoarseData
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import get_coarse_structure, load_coarse_structure
from packs import defnames
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights

import numpy as np
from typing import Sequence
from shapely import geometry
from scipy.sparse.linalg import spsolve


def get_perm_diag(value):
    return np.array([[value, 0], [0, value]])

def set_permeability_brazil(fine_mesh_path, fine_properties: MeshProperty, typek='barrier', export_permfield=True, update_permfield=True) -> None:
    typeks = ['barrier', 'channel']
    assert typek in typeks

    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return

    k1 = 1.0
    k2 = 1e-3
    k3 = 1e3

    faces_in_losangle = fine_properties['faces_in_losangle']
    faces = fine_properties['faces']

    permeability = np.zeros((faces.shape[0], 2, 2))
    permeability[:] = get_perm_diag(k1)
    if typek == typeks[0]:
        permeability[faces_in_losangle] = permeability[faces_in_losangle]*k2
    else:
        permeability[faces_in_losangle] = permeability[faces_in_losangle]*k3
    
    fine_properties.insert_or_update_data({
        tag_preprocess: permeability
    })

    if export_permfield:
        mesh_data = MeshData(mesh_path=fine_mesh_path)
        mesh_data.create_tag('permeability')
        mesh_data.insert_tag_data(
            'permeability',
            fine_properties['permeability'][:, 1, 1],
            elements_type='faces'
        )
        mesh_data.export_all_elements_type_to_vtk('permfield', element_type='faces')

def define_faces_in_losangle(fine_properties: MeshProperty) -> None:
    tag_preprocess = 'faces_in_losangle'
    if fine_properties.verify_name_in_data_names(tag_preprocess):
        return
    
    Lx = 1.5
    Ly = 1

    d1 = 1.2
    d2 = 0.75

    x1 = (Lx-d1)/2
    x2 = x1 + d1/2
    x3 = x1 + d1
    x4 = x2

    y1 = Ly/2
    y2 = (Ly-d2)/2
    y3 = y1
    y4 = y2 + d2

    losangle = geometry.Polygon([
        (x1, y1),
        (x2, y2),
        (x3, y3),
        (x4, y4)
    ])

    poly = losangle
    faces_centroids = fine_properties['faces_centroids']
    points_list = geometry.MultiPoint(faces_centroids)
    test = np.array([poly.contains(i) for i in points_list.geoms])
    
    faces_in_losangle = fine_properties['faces'][test]

    fine_properties.insert_or_update_data({
        tag_preprocess: faces_in_losangle
    })

def define_coarse_structure(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation, level=1, update=True) -> Sequence[PrimalCoarseData]:
    fp = fine_mesh_properties
    if update is True:
        coarse_struct = get_coarse_structure(
            1,
            fp[defnames.get_primal_id_name_by_level(level)],
            fp['faces'],
            fp['adjacencies'],
            fp['edges'],
            fp['nodes_of_edges'],
            fp['bool_boundary_edges'],
            fp['bool_boundary_nodes'],
            fp['nodes'],
            fp['nodes_weights'],
            fp['nodes_of_nodes'],
            fp['edges_of_nodes'],
            fp['faces_of_nodes'],
            fp['nodes_centroids'],
            fp['faces_centroids'],
            fp['permeability'],
            fp['unitary_normal_edges'],
            fp[defnames.get_dual_id_name_by_level(1)],
            fp.edges_dim,
            lsds
        )
    else:
        coarse_struct = load_coarse_structure(
            level,
            fp[defnames.get_primal_id_name_by_level(level)]
        )
    
    return coarse_struct

def set_weights_nodes_cstruct(cstruct: PrimalCoarseData):

    fine_properties = cstruct
    weights = get_gls_nodes_weights(**fine_properties.get_all_data())
    return weights

def update_xi_params(xi_params, total_mobility_edges):
    xi_params_new = xi_params.copy()
    xi_params_new[:] = xi_params*total_mobility_edges[:, np.newaxis]
    return xi_params_new

def set_fine_transmissibility_biphasic(fine_mesh_properties: MeshProperty, bc: BoundaryConditions, lsds: LsdsFluxCalculation):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )
    return resp

def update_fine_flux(
        coarse_struct: Sequence[PrimalCoarseData],
        edges_flux: np.ndarray,
        ms_pressure: np.ndarray,
        total_mobility_edges: np.ndarray,
        lsds: LsdsFluxCalculation,
        nodes_pressure: np.ndarray
    ):

    for cstruct in coarse_struct:
        global_edges = cstruct['map_edges']
        bool_boundary_edges = cstruct['bool_boundary_edges']
        local_flux_presc = edges_flux[global_edges[bool_boundary_edges]]
        neumann_edges = cstruct['edges'][bool_boundary_edges]
        bc = BoundaryConditions()
        bc.set_boundary('neumann_edges', neumann_edges, local_flux_presc)
        # bc.set_boundary('neumann_edges', np.array([]), np.array([]))

        local_vertice = cstruct['faces'][cstruct['dual_id']==defnames.dual_ids('vertice_id')]
        bc.set_boundary('dirichlet_volumes', local_vertice, ms_pressure[local_vertice])
        # bc.set_boundary('dirichlet_volumes', np.array([]), np.array([]))
        bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))

        # bool_boundary_nodes = cstruct['bool_boundary_nodes']
        # mapbnodes = cstruct['map_nodes'][bool_boundary_nodes]
        # bnodes = cstruct['nodes'][bool_boundary_nodes]
        # bc.set_boundary('dirichlet_nodes', bnodes, nodes_pressure[mapbnodes])

        edges_multiplier = total_mobility_edges[global_edges].copy()
        ########################
        # edges_multiplier[:] = 1
        #########################
        cstruct.update_data({
            'neumann_edges': neumann_edges,
            'neumann_edges_value': local_flux_presc,
            'edges_multiplier': edges_multiplier,
        })
        
        boundary_nodes_weights = set_weights_nodes_cstruct(cstruct)
        local_nodes_weights = cstruct['nodes_weights_internal'].copy()
        local_nodes_weights = np.hstack([local_nodes_weights, boundary_nodes_weights['nodes_weights']])

        cstruct.insert_or_update_data({
            'xi_params': update_xi_params(cstruct['xi_params_backup'], total_mobility_edges[global_edges]),
            'nodes_weights': local_nodes_weights,
            'neumann_weights': boundary_nodes_weights['neumann_weights']
        })

        lt = set_fine_transmissibility_biphasic(cstruct, bc, lsds)
        local_pressure = spsolve(lt['transmissibility'], lt['source'])
        
        local_edges_flux = lsds.get_edges_flux(
            bc,
            local_pressure,
            cstruct['xi_params'],
            local_nodes_weights,
            cstruct['nodes_of_edges'],
            cstruct['adjacencies'],
            cstruct['neumann_weights']
        )

        # if cstruct['coarse_id'][0] == 1:
        #     import pdb; pdb.set_trace()

        cstruct.insert_or_update_data({
            'edges_flux': local_edges_flux
        })

        bool_internal_edges = ~bool_boundary_edges
        # edges_flux[cstruct['map_edges'][bool_internal_edges]] = local_edges_flux[bool_internal_edges]
        edges_flux[cstruct['map_edges']] = local_edges_flux