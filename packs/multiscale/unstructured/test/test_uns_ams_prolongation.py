from .create_primal_test import get_fine_mesh_path_and_mesh_properties_name_for_test, get_coarse_mesh_path_and_mesh_properties_name_for_test
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio_if_not_exists, _create_flying_mesh
from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d import create_dual
# from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d_v2 import create_dual
from packs import defnames
from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.multiscale.unstructured.operators.prolongation.ams import Unstructured2DAmsOperator
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.utils.profile_functions import profile
from packs.utils.multiscale_methods import print_adm_interfaces_2d
from packs.adm.adm_operators import Adm
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import create_dual_interaction_regions
from packs.multiscale.unstructured.operators.prolongation.get_op_from_amsu import update_global_op_from_amsu

import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp
import os
from typing import Sequence


def create_primal_ids(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty):
    key1 = defnames.get_primal_id_name_by_level(1)
    if not fine_mesh_properties.verify_name_in_data_names(key1):
    
        fine_primal_ids = create_coarse_volumes(
            faces_id_level0=fine_mesh_properties['faces'],
            faces_centroids_level0=fine_mesh_properties['faces_centroids'],
            faces_ids_level1=coarse_mesh_properties['faces'],
            nodes_centroids_level1=coarse_mesh_properties['nodes_centroids'],
            nodes_of_faces_level1=coarse_mesh_properties['nodes_of_faces'],
            adjacencies_level0=fine_mesh_properties['adjacencies'],
            faces_of_faces_level0=fine_mesh_properties.faces_of_faces,
            faces_centroids_level1=coarse_mesh_properties['faces_centroids'],
            faces_of_faces_level1=coarse_mesh_properties.faces_of_faces,
            level=1
        )

        fine_mesh_properties.insert_or_update_data(
            fine_primal_ids
        )

        fine_mesh_properties.export_data()

def create_dual_ids(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty):
    key1 = defnames.get_dual_id_name_by_level(1)
    if key1 not in fine_mesh_properties.keys():

        dual_data = create_dual(fine_mesh_properties=fine_mesh_properties, coarse_mesh_properties=coarse_mesh_properties, level=1)

        fine_mesh_properties.insert_or_update_data(
            dual_data
        )

        fine_mesh_properties.export_data()

def define_boundary_conditions(fine_mesh_properties: MeshProperty):
    bc_name = 'test1_ams'
    bc = BoundaryConditions(name=bc_name)

    dirichlet_edges_1 = fine_mesh_properties['Inflow']
    dirichlet_edges_0 = fine_mesh_properties['Outflow']
    walls_edges = fine_mesh_properties['Walls']

    dirichlet_nodes1 = np.unique(
        fine_mesh_properties['nodes_of_edges'][
            dirichlet_edges_1
        ])
    
    dirichlet_nodes0 = np.unique(
        fine_mesh_properties['nodes_of_edges'][
            dirichlet_edges_0
        ])
    
    bc_nodes = np.concatenate([dirichlet_nodes1, dirichlet_nodes0])
    nodes_values = np.concatenate([
        np.repeat(10.0, dirichlet_nodes1.shape[0]),
        np.repeat(1.0, dirichlet_nodes0.shape[0])
    ])

    bc.set_boundary('dirichlet_nodes', bc_nodes, nodes_values)

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    return bc

def update_permeability(fine_mesh_properties: MeshProperty):
    np.random.seed(5)
    n_faces = fine_mesh_properties['faces'].shape[0]
    permeability = np.zeros((n_faces, 2, 2))
    v1 = np.random.randint(0, 6, size=n_faces)
    perms = np.float_power(10, -v1)
    permeability[:, 0, 0] = perms
    permeability[:, 1, 1] = perms

    fine_mesh_properties.insert_or_update_data({
        'permeability': permeability
    })

def func1(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation):
    transm = lsds.mount_transmissibility_matrix_without_bc(**fine_mesh_properties.get_all_data())
    return transm

def func2(lsds: LsdsFluxCalculation, bc: BoundaryConditions, fine_mesh_properties: MeshProperty):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )
    return resp

# @profile
def func3(ams_prolongation: Unstructured2DAmsOperator, transmissibility_without_bc, diagonal_term):
    local_matrices = ams_prolongation.get_local_transmissibility_matrix(
        ams_prolongation['dual_volumes'],
        transmissibility_without_bc,
        diagonal_term
    )
    return local_matrices




def export_op(mesh_path, OP_AMS):

    # flying_mesh_path = _create_flying_mesh(mesh_path)
    flying_mesh_path = mesh_path


    mesh_data = MeshData(mesh_path=flying_mesh_path)
    all_data = sp.find(OP_AMS)
    lines = all_data[0]
    cols = all_data[1]
    data = all_data[2]

    elements = []
    data_array = []

    cids = np.unique(cols)

    for cid in cids:
        test = cols == cid
        elements.append(lines[test])
        data_array.append(data[test])

    mesh_data.insert_array_tag_data(
        'OP',
        data_array,
        'faces',
        elements
    )

    mesh_data.export_all_elements_type_to_vtk(
        'op_amsU',
        'faces'
    )
    
def export_adm_levels(mesh_path, fine_levels):

    flying_mesh_path = _create_flying_mesh(mesh_path)
    mesh_data = MeshData(mesh_path=flying_mesh_path)
    mesh_data.create_tag('adm_levels', data_type='int')
    mesh_data.insert_tag_data('adm_levels', fine_levels, 'faces')
    mesh_data.export_all_elements_type_to_vtk('adm_test', 'faces')



def run():
    lsds = LsdsFluxCalculation()

    fine_mesh_path, fine_mesh_properties_name, fine_mesh_path_v4 = get_fine_mesh_path_and_mesh_properties_name_for_test()
    coarse_mesh_path, coarse_mesh_properties_name = get_coarse_mesh_path_and_mesh_properties_name_for_test()

    # fine_mesh_properties = create_meshproperties_from_meshio_if_not_exists(fine_mesh_path, fine_mesh_properties_name)
    # coarse_mesh_properties = create_meshproperties_from_meshio_if_not_exists(coarse_mesh_path, coarse_mesh_properties_name)

    fine_mesh_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)
    coarse_mesh_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    # mesh_data = MeshData(mesh_path = fine_mesh_path)
    # mesh_data.create_tag('permeability')
    # perms = fine_mesh_properties.permeability[:, 0, 0]
    # mesh_data.insert_tag_data('permeability', perms, elements_type='faces')
    # mesh_data.export_all_elements_type_to_vtk('perm_field', element_type='faces')


    create_primal_ids(fine_mesh_properties, coarse_mesh_properties)
    create_dual_ids(fine_mesh_properties, coarse_mesh_properties)

    # update_permeability(fine_mesh_properties)
    
    ams_prolongation = Unstructured2DAmsOperator()

    ams_prolongation.insert_ams_data(
        {
            # defnames.dual_volumes_str: fine_mesh_properties[defnames.get_dual_volumes_name_by_level(1)],
            defnames.dual_volumes_str: np.array([fine_mesh_properties['faces']]),
            defnames.fine_primal_id: fine_mesh_properties[defnames.get_primal_id_name_by_level(1)],
            defnames.fine_dual_id: fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
        }
    )

    ams_prolongation.preprocess_ams_data()

    bc = define_boundary_conditions(fine_mesh_properties)

    fine_mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    if not fine_mesh_properties.verify_name_in_data_names('nodes_weights'):
        fine_mesh_properties.insert_or_update_data(
            {'nodes_to_calculate': fine_mesh_properties['nodes']}
        )
        fine_mesh_properties.insert_or_update_data(
            get_gls_nodes_weights(**fine_mesh_properties.get_all_data())
        )
        fine_mesh_properties.export_data()
    
    if not fine_mesh_properties.verify_name_in_data_names('xi_params'):
        fine_mesh_properties.insert_or_update_data(
            lsds.get_all_edges_flux_params(**fine_mesh_properties.get_all_data())
        )
        fine_mesh_properties.export_data()
    
    # transm = lsds.mount_transmissibility_matrix_without_bc(**fine_mesh_properties.get_all_data())
    transm = func1(fine_mesh_properties, lsds)
    monotone_transm = ams_prolongation.get_monotone_matrix(
        transm['transmissibility_without_bc'],
        epsilon=0.0001,
        w=0.97
    )

    
    # resp = lsds.mount_transmissibility_matrix(
    #     bc,
    #     **fine_mesh_properties.get_all_data()
    # )

    resp = func2(lsds, bc, fine_mesh_properties)

    level_str = defnames.level_str(1)

    interaction_regions = create_dual_interaction_regions(
        fine_mesh_properties[defnames.get_dual_interation_region_name_by_level(1)],
        fine_mesh_properties[defnames.vertices_selected + level_str],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)][fine_mesh_properties[defnames.vertices_selected + level_str]],
        fine_mesh_properties[defnames.boundary_dual_interaction + level_str],
        fine_mesh_properties[defnames.internal_dual_path + level_str],
        fine_mesh_properties[defnames.dual_initial_ccs + level_str],
        global_transmissibility=transm['transmissibility_without_bc'],
        global_diagonal_term=np.zeros(resp['source'].shape[0])
    )  

    OP_AMS = ams_prolongation.get_global_op(coarse_mesh_properties['faces'], fine_mesh_properties['faces'])

    # ams_prolongation.insert_data_in_global_OP(
    #     ams_prolongation[ams_prolongation.dual_volumes_str],
    #     monotone_transm,
    #     np.zeros(resp['source'].shape[0]),
    #     OP_AMS
    # )

    update_global_op_from_amsu(interaction_regions, OP_AMS)

    OR_AMS = ams_prolongation.get_finite_volume_restriction_operator(
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)]
    )

    export_op(fine_mesh_path, OP_AMS)
    import pdb; pdb.set_trace()

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

    ### adm
    adm = Adm()
    fine_levels = np.full(len(fine_mesh_properties['faces']), -1)
    boundary_coarse_faces = np.unique(np.concatenate(coarse_mesh_properties.faces_of_nodes[
        coarse_mesh_properties['bool_boundary_nodes']
        ]))
    fine_levels[
        np.isin(fine_mesh_properties[defnames.get_primal_id_name_by_level(1)], boundary_coarse_faces)
    ] = 0

    list_primal_ids = [
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)]
    ]

    list_dual_ids = [
        fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    ]

    adm.update_levels(
        list_primal_ids,
        fine_levels,
        fine_mesh_properties['faces'],
        2,
        fine_mesh_properties.faces_of_faces_by_nodes
    )

    # export_adm_levels(fine_mesh_path, fine_levels)
    # flying_mesh_path = _create_flying_mesh(fine_mesh_path)
    # print_adm_interfaces_2d(
    #     fine_mesh_properties,
    #     flying_mesh_path,
    #     fine_levels,
    #     'adm_edges'
    # )

    OP_adm, OR_adm = adm.get_adm_prolongation_operator(
        [OP_AMS],
        [OR_AMS],
        fine_levels,
        list_primal_ids,
        list_dual_ids,
        fine_mesh_properties['faces'],
        2
    )

    T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    Q_adm = OR_adm*resp['source']
    P_adm = spsolve(T_adm.tocsc(), Q_adm)
    P_prol = OP_adm*P_adm

    error = np.absolute(pressure - P_prol)
    relative_error = (error/pressure)*100

    mesh_data = MeshData(mesh_path = os.path.join('remove', 'flying_mesh.h5m'))
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, elements_type='faces')

    mesh_data.create_tag('adm_prol_pressure')
    mesh_data.insert_tag_data('adm_prol_pressure', P_prol, elements_type='faces')

    mesh_data.create_tag('absolute_error')
    mesh_data.insert_tag_data('absolute_error', error, elements_type='faces')

    mesh_data.create_tag('relative_error')
    mesh_data.insert_tag_data('relative_error', relative_error, elements_type='faces')

    mesh_data.export_all_elements_type_to_vtk('adm_solution_classic_global', element_type='faces')






    import pdb; pdb.set_trace()

    # edges_flux = lsds.get_edges_flux(
    #     bc,
    #     pressure,
    #     fine_mesh_properties.xi_params,
    #     fine_mesh_properties.nodes_weights,
    #     fine_mesh_properties.nodes_of_edges,
    #     fine_mesh_properties.adjacencies,
    #     fine_mesh_properties['neumann_weights']      
    # )

    # faces_flux = lsds.get_faces_flux(
    #     edges_flux,
    #     fine_mesh_properties.adjacencies,
    #     fine_mesh_properties.bool_boundary_edges
    # )

    # import pdb; pdb.set_trace()


    







    # import pdb; pdb.set_trace()




    return ams_prolongation, fine_mesh_properties, coarse_mesh_properties




    # key_str = defnames.get_dual_id_name_by_level(level=1)
    # data = dual_data.get(key_str)

    # flying_fine_mesh_path = _create_flying_mesh(fine_mesh_path)
    # mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    # mesh_data.create_tag(key_str, data_type='int')
    # mesh_data.insert_tag_data(key_str, data, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    # mesh_data.export_only_the_elements(key_str, element_type='faces', elements_array=fine_mesh_properties['faces'])

    # dual_volumes_name = defnames.get_dual_volumes_name_by_level(1)
    # dual_volumes = dual_data[dual_volumes_name]

    # mesh_data.export_list_elements_array_data(dual_volumes_name, 'faces', dual_volumes)

    # interaction_regions_name = defnames.get_dual_interation_region_name_by_level(1)
    # regions = dual_data[interaction_regions_name]

    # mesh_data.export_list_elements_array_data(interaction_regions_name, 'faces', regions)




    # flying_coarse_mesh_path = _create_flying_mesh(coarse_mesh_path)
    # mesh_data = MeshData(mesh_path=flying_coarse_mesh_path)   
    # mesh_data.create_tag('id', data_type='int')
    # mesh_data.insert_tag_data('id', coarse_mesh_properties['faces'], elements_type='faces', elements_array=coarse_mesh_properties['faces'])
    # mesh_data.export_only_the_elements('coarse_ids_test', element_type='faces', elements_array=coarse_mesh_properties['faces'])


    print('fim')



