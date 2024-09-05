from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d import create_dual
from packs import defnames
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.utils import utils_old
from packs.multiscale.unstructured.operators.prolongation.ams import Unstructured2DAmsOperator
from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import create_dual_interaction_regions
from packs.multiscale.unstructured.operators.prolongation.get_op_from_amsu import update_global_op_from_amsu
from packs.adm.non_uniform import fine_level_from_alpha, fine_level_from_beta
from packs.fim_nu_adm.packs.processor import nu_adm_funcs
from packs.utils.multiscale_methods import print_adm_interfaces_2d
from packs.multiscale.unstructured.test.test_uns_ams_prolongation import export_adm_levels
from packs.utils.multiscale_methods import print_fine_interfaces_coarse_mesh_2d


import os
from shapely import geometry
from typing import List, Union
import numpy as np
from scipy.sparse.linalg import spsolve, gmres, cg, bicgstab, spilu, LinearOperator
import scipy.sparse as sp
import pandas as pd

def get_perm_diag(value):
    return np.array([[value, 0], [0, value]])

def get_properties():
    rel_path = os.path.join(
        defpaths.unstructured_coarse_test_mesh_folder,
        'brazil'
    )
    # fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
    # fine_mesh_properties_name = 'brazilf' 
    # fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    fine_mesh_path = os.path.join(rel_path, 'brazilfquad.msh')
    fine_mesh_properties_name = 'brazilfquad' 
    fine_mesh_path_v4 = os.path.join(rel_path, 'brazilfquad_v4.msh')

    # fine_mesh_path = os.path.join(rel_path, 'brazilf_test.msh')
    # fine_mesh_properties_name = 'brazilf_test' 
    # fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_test_v4.msh')

    # coarse_mesh_path = os.path.join(rel_path, 'brazilC.msh')
    # coarse_mesh_properties_name = 'brazilC1'

    coarse_mesh_path = os.path.join(rel_path, 'brazilC2.msh')
    coarse_mesh_properties_name = 'brazilC2'

    # coarse_mesh_path = os.path.join(rel_path, 'brazilC3.msh')
    # coarse_mesh_properties_name = 'brazilC3'

    # coarse_mesh_path = os.path.join(rel_path, 'brazilC4.msh')
    # coarse_mesh_properties_name = 'brazilC4'

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)
    coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    return fine_properties, coarse_properties, fine_mesh_path, coarse_mesh_path

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

def define_boundary(fine_properties: MeshProperty) -> None:

    data_to_update = dict()

    nodes = fine_properties['nodes']
    nodes_centroids = fine_properties['nodes_centroids']
    edges = fine_properties['edges']
    edges_centroids = fine_properties.edges_centroids

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    center = np.array([xmax/2, ymax/2])
    R = 0.2

    delta = fine_properties.edges_dim.min()/10

    dists = np.linalg.norm(nodes_centroids - center, axis=1)
    nodes_internal_boundary = nodes[
        dists <= R + delta
    ]

    data_to_update.update({'nodes_internal_boundary': nodes_internal_boundary})

    

    inflow = edges[edges_centroids[:, 0] < xmin + delta]
    outflow = edges[edges_centroids[:, 0] > xmax - delta]

    walls = np.concatenate([
        edges[edges_centroids[:, 1] < ymin + delta],
        edges[edges_centroids[:, 1] > ymax - delta]
    ])

    data_to_update.update({
        'Inflow': inflow,
        'Outflow': outflow,
        'Walls': walls
    })

    fine_properties.insert_or_update_data(data_to_update)




def set_permeability(fine_mesh_path, fine_properties: MeshProperty, typek='barrier', export_permfield=True, update_permfield=True) -> None:
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

def set_boundary_conditions(fine_properties: MeshProperty) -> BoundaryConditions:
    
    bc = BoundaryConditions()

    k = 1

    nodes = fine_properties['nodes']
    nodes_centroids = fine_properties['nodes_centroids']
    edges = fine_properties['edges']
    edges_centroids = fine_properties.edges_centroids

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    R = 0.2*k

    delta = fine_properties.edges_dim.min()/10

    center = np.array([xmax/2, ymax/2])

    nodes_p1 = nodes[nodes_centroids[:,0] < xmin + delta]
    nodes_p0 = nodes[nodes_centroids[:,0] > xmax - delta]

    dist_nodes = np.linalg.norm(
        nodes_centroids - center,
        axis=1
    )

    nodes_p05 = nodes[dist_nodes < R + delta]

    edges_ymax = edges[edges_centroids[:, 1] > ymax - delta]
    edges_ymin = edges[edges_centroids[:, 1] < ymin + delta]

    bc_nodes = np.concatenate([nodes_p1, nodes_p0, nodes_p05])
    nodes_values = np.concatenate([
        np.repeat(1.0, nodes_p1.shape[0]),
        np.repeat(0.0, nodes_p0.shape[0]),
        np.repeat(0.5, nodes_p05.shape[0])
    ])

    bc.set_boundary('dirichlet_nodes', bc_nodes, nodes_values)

    walls_edges = np.unique(np.concatenate([edges_ymax, edges_ymin]))

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    return bc

def set_weights_nodes(fine_properties: MeshProperty, update=True):

    calculate = False

    if update is True:
        calculate = True
        
    if not fine_properties.verify_name_in_data_names('nodes_weights'):
        calculate=True

    if calculate is True:
        fine_properties.insert_or_update_data(
            {'nodes_to_calculate': fine_properties['nodes']}
        )
        weights = get_gls_nodes_weights(**fine_properties.get_all_data())
        fine_properties.insert_or_update_data(weights)

        lsds = LsdsFluxCalculation()
        fine_properties.insert_or_update_data(
            lsds.get_all_edges_flux_params(**fine_properties.get_all_data())
        )

        fine_properties.export_data()

def create_primal_ids(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty, update=True):
    create = False
    
    key1 = defnames.get_primal_id_name_by_level(1)

    if update is True:
        create = True

    if not fine_mesh_properties.verify_name_in_data_names(key1):
        create = True
    
    if create is True:
    
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
            level=1,
            edges_ids_level0=fine_mesh_properties['edges'],
            bool_boundary_edges_level0=fine_mesh_properties['bool_boundary_edges'],
            edges_centroids_level0=fine_mesh_properties.edges_centroids
        )

        fine_mesh_properties.insert_or_update_data(
            fine_primal_ids
        )

        fine_mesh_properties.export_data()

def export_primal_ids(fine_mesh_path, fine_mesh_properties: MeshProperty, coarse_mesh_path, export=True):
    if export is True:
        pass
    else:
        return
    
    key1 = defnames.get_primal_id_name_by_level(1)
    primal_id = fine_mesh_properties[key1]

    key_str = key1
    data = primal_id

    flying_fine_mesh_path = fine_mesh_path
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(key_str, data_type='int')
    mesh_data.insert_tag_data(key_str, data, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements(key_str, element_type='faces', elements_array=fine_mesh_properties['faces'])
    

    print_fine_interfaces_coarse_mesh_2d(
        fine_mesh_properties,
        flying_fine_mesh_path,
        1,
        'edges_selected_2'
    )

    coarse_mesh_data = MeshData(mesh_path=coarse_mesh_path)
    coarse_mesh_data.export_all_elements_type_to_vtk('background_coarse_mesh', element_type='faces')

def create_dual_ids(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty, update=True, dual_type=1):
    create = False
    key1 = defnames.get_dual_id_name_by_level(1)

    if update is True:
        create = True
    if not fine_mesh_properties.verify_name_in_data_names(key1):
        create = True
    
    if create is True:

        dual_data = create_dual(fine_mesh_properties=fine_mesh_properties, coarse_mesh_properties=coarse_mesh_properties, level=1, dual_type=dual_type)

        fine_mesh_properties.insert_or_update_data(
            dual_data
        )

        fine_mesh_properties.export_data()

def export_dual_ids(fine_mesh_path, fine_mesh_properties: MeshProperty, export=True):
    if export is True:
        pass
    else:
        return

    key_str = defnames.get_dual_id_name_by_level(level=1)
    data = fine_mesh_properties[key_str]

    flying_fine_mesh_path = fine_mesh_path
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(key_str, data_type='int')
    mesh_data.insert_tag_data(key_str, data, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements(key_str, element_type='faces', elements_array=fine_mesh_properties['faces'])

    dual_volumes_name = defnames.get_dual_volumes_name_by_level(1)
    dual_volumes = fine_mesh_properties[dual_volumes_name]

    mesh_data.export_list_elements_array_data(dual_volumes_name, 'faces', dual_volumes)

    interaction_regions_name = defnames.get_dual_interation_region_name_by_level(1)
    regions = fine_mesh_properties[interaction_regions_name]

    mesh_data.export_list_elements_array_data(interaction_regions_name, 'faces', regions)

    dual_boundarys_name = 'dual_boundary'
    boundarys = fine_mesh_properties[defnames.boundary_dual_interaction + defnames.level_str(1)]

    mesh_data.export_list_elements_array_data(dual_boundarys_name, 'faces', boundarys)

def set_fine_transmissibility_without_bc(save_fine_transm_without_bc, matrix_path, fine_mesh_properties: MeshProperty, fine_transm_without_bc_name):
    transm = dict()
    
    if save_fine_transm_without_bc is True:
        lsds = LsdsFluxCalculation()
        transm.update(func1(fine_mesh_properties, lsds))
        utils_old.save_matrix(matrix_path, fine_transm_without_bc_name, matrix=transm.get('transmissibility_without_bc'))
    else:
        transm.update({
            'transmissibility_without_bc': utils_old.load_matrix(
                matrix_path,
                fine_transm_without_bc_name
            )
        })
    
    return transm

def func1(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation):
    transm = lsds.mount_transmissibility_matrix_without_bc(**fine_mesh_properties.get_all_data())
    return transm

def set_monotone_transm(save_monotone_transm, transm: dict, w, matrix_path, monotone_transm_name):
    ams_prolongation = Unstructured2DAmsOperator()
    if save_monotone_transm is True:
        monotone_transm = ams_prolongation.get_monotone_matrix(
            transm['transmissibility_without_bc'],
            epsilon=0.001,
            w=w
        )
        
        # enhanced = Enhanced()
        # monotone_transm = enhanced.get_enhanced_matrix(transm['transmissibility_without_bc'])
        
        utils_old.save_matrix(matrix_path, monotone_transm_name, matrix=monotone_transm)
    else:
        monotone_transm = utils_old.load_matrix(
            matrix_path,
            monotone_transm_name
        )
    
    return monotone_transm

def func2(lsds: LsdsFluxCalculation, bc: BoundaryConditions, fine_mesh_properties: MeshProperty):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )
    return resp

def set_fine_transmissibility(save_fine_transmissibility, fine_mesh_properties: MeshProperty, bc: BoundaryConditions, matrix_path, fine_transmissibility_name):
    resp = dict()
    lsds = LsdsFluxCalculation()

    if save_fine_transmissibility is True:

        resp.update(func2(lsds, bc, fine_mesh_properties))
        utils_old.save_matrix(matrix_path, fine_transmissibility_name, matrix=resp['transmissibility'])
        fine_mesh_properties.insert_or_update_data({
            'fine_source': resp['source']
        })
        fine_mesh_properties.export_data()
    else:

        resp.update(
            {
                'transmissibility': utils_old.load_matrix(
                    matrix_path,
                    fine_transmissibility_name
                ),
                'source': fine_mesh_properties['fine_source']
            }
        )
    
    return resp

def get_OR_AMS(fine_mesh_properties: MeshProperty):
    ams_prolongation = Unstructured2DAmsOperator()
    OR_AMS = ams_prolongation.get_finite_volume_restriction_operator(
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)]
    )
    return OR_AMS

def get_op( 
        save_op, 
        fine_mesh_properties: MeshProperty, 
        coarse_mesh_properties: MeshProperty,
        op_name,
        matrix_path,
        op_toget,
        monotone_transm,
        resp,
        level_str
    ):

    ams_prolongation = Unstructured2DAmsOperator()
    if save_op is True:
        OP_AMS = ams_prolongation.get_global_op(coarse_mesh_properties['faces'], fine_mesh_properties['faces'])

        if op_toget == 'AMS-U':

                interaction_regions = create_dual_interaction_regions(
                    fine_mesh_properties[defnames.get_dual_interation_region_name_by_level(1)],
                    fine_mesh_properties[defnames.vertices_selected + level_str],
                    fine_mesh_properties[defnames.get_primal_id_name_by_level(1)][fine_mesh_properties[defnames.vertices_selected + level_str]],
                    fine_mesh_properties[defnames.internal_dual_path + level_str],
                    fine_mesh_properties[defnames.boundary_dual_interaction + level_str],
                    fine_mesh_properties[defnames.dual_initial_ccs + level_str],
                    global_transmissibility=monotone_transm,
                    global_diagonal_term=np.zeros(resp['source'].shape[0]),
                    dual_id=fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
                )

                OP_AMS = update_global_op_from_amsu(interaction_regions, OP_AMS)
        
        utils_old.save_matrix(matrix_path, op_name, OP_AMS)
        # export_op(fine_mesh_path, OP_AMS, op_name)
    else:
       OP_AMS = utils_old.load_matrix(matrix_path, op_name)
    
    return OP_AMS

def define_new_fine_levels_v0(
    fine_mesh_properties: MeshProperty,
    bc: BoundaryConditions
) -> np.ndarray:
    """
    Nenhum volume na malha fina
    """
    return np.array([])

def define_new_fine_levels_v1(
        fine_mesh_properties: MeshProperty,
        bc: BoundaryConditions
) -> np.ndarray:
    
    """
    Apenas as duais com prescricao na malha fina usando a dual tipo 1
    """

    faces_of_nodes = fine_mesh_properties['faces_of_nodes']
    dual_volumes = fine_mesh_properties['dual_volumes_level1']
    dual_id = fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    
    nodes_pressure_presc = bc['dirichlet_nodes']['id']
    faces_of_nodes_presc = np.unique(
        np.concatenate(faces_of_nodes[nodes_pressure_presc])
    )
    faces_of_nodes_presc = faces_of_nodes_presc[dual_id[faces_of_nodes_presc] == defnames.dual_ids('face_id')]
    
    boundary_faces = faces_of_nodes_presc
    dual_in_boundary = []
    for dual in dual_volumes:
        if np.any(np.isin(dual, boundary_faces)):
            dual_in_boundary.append(dual)
    
    dual_in_boundary = np.unique(np.concatenate(dual_in_boundary))
    return dual_in_boundary

def define_new_fine_levels_v2(
        fine_mesh_properties: MeshProperty,
        bc: BoundaryConditions
) -> np.ndarray:
    
    """
    Apenas as faces com prescricao na malha fina
    """

    faces_of_nodes = fine_mesh_properties['faces_of_nodes']
    dual_volumes = fine_mesh_properties['dual_volumes_level1']
    dual_id = fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    
    nodes_pressure_presc = bc['dirichlet_nodes']['id']
    faces_of_nodes_presc = np.unique(
        np.concatenate(faces_of_nodes[nodes_pressure_presc])
    )
    # faces_of_nodes_presc = faces_of_nodes_presc[dual_id[faces_of_nodes_presc] == defnames.dual_ids('face_id')]
    
    # boundary_faces = faces_of_nodes_presc
    # dual_in_boundary = []
    # for dual in dual_volumes:
    #     if np.any(np.isin(dual, boundary_faces)):
    #         dual_in_boundary.append(dual)
    
    # dual_in_boundary = np.unique(np.concatenate(dual_in_boundary))
    return faces_of_nodes_presc

def write_results(
        l2_error, 
        l2_relative_error, 
        max_abs_error, 
        max_abs_relative_error,
        dual_type,
        fine_mesh_name,
        coarse_mesh_name,
        number_fine_vols,
        number_coarse_vols,
        perm_type,
        fine_levels_setup,
    ):

    new_data = np.array([
        l2_error,
        l2_relative_error,
        max_abs_error,
        max_abs_relative_error,
        dual_type,
        fine_mesh_name,
        coarse_mesh_name,
        number_fine_vols,
        number_coarse_vols,
        perm_type
    ], dtype='O')

    file_path = os.path.join(defpaths.flying, 'sim_results.csv')

    data_types = [np.float64, np.float64, np.float64, np.float64, np.int, 'O', 'O', np.int, np.int, 'O', np.int]

    header = np.array([
            'l2_error',
            'l2_relative_error',
            'max_abs_error',
            'max_abs_relative_error',
            'dual_type',
            'fine_mesh_name',
            'coarse_mesh_name',
            'number_fine_vols',
            'number_coarse_vols',
            'perm_type',
            'fine_levels_setup'
        ], dtype='<U30')

    if os.path.exists(file_path):
        pass
    else:
        first_line = np.array([np.array([0.0], dtype=t) for t in data_types], dtype='O').T
        df = pd.DataFrame(first_line, columns=header)
        for i in range(header.shape[0]):
            df[header[i]] = df[header[i]].astype(data_types[i])
        df.to_csv(file_path, index=False)
    
    my_types = dict()
    for i in range(header.shape[0]):
        my_types.update({header[i]: data_types[i]})
        
    df = pd.read_csv(file_path, dtype=my_types)

    new_df = {
        'l2_error': l2_error,
        'l2_relative_error': l2_relative_error,
        'max_abs_error': max_abs_error,
        'max_abs_relative_error': max_abs_relative_error,
        'dual_type': dual_type,
        'fine_mesh_name': fine_mesh_name,
        'coarse_mesh_name': coarse_mesh_name,
        'number_fine_vols': number_fine_vols,
        'number_coarse_vols': number_coarse_vols,
        'perm_type': perm_type,
        'fine_levels_setup': fine_levels_setup
    }
    df.loc[len(df)] = new_df
    df.reset_index(drop=True, inplace=True)
    subset_to_drop = header[4:]
    df.drop_duplicates(subset=subset_to_drop, inplace=True)
    df.to_csv(file_path, index=False)

    





def run4():
    matrix_path = 'matrices.h5'
    fine_transm_without_bc_name = 'fine_transm_without_bc'
    save_fine_transm_without_bc = False
    save_monotone_transm = False
    save_fine_transmissibility = False
    save_op = False
    w = 0
    monotone_transm_name = 'monotone_transm_w_0'
    fine_transmissibility_name = 'fine_transmissibility'
    op_toget = 'AMS-U'
    op_name = 'AMS_U_w_0_dual2'
    level_str = defnames.level_str(1)
    alpha_lim_finescale = 0.5
    beta_lim = 3
    export_adm_levels_file = True
    bool_export_primal_id = False
    bool_export_dual_id = False
    my_dual_type = 1
    perm_type = 'barrier'
    update_nodes_weights = False
    export_permfield = False
    update_permfield = False
    fine_level_setup = 2

    fp, cp, fine_mesh_path, coarse_mesh_path = get_properties()
    mesh_data = MeshData(mesh_path=coarse_mesh_path)
    mesh_data.export_all_elements_type_to_vtk('background_coarse_mesh', 'faces')
    define_faces_in_losangle(fp)
    set_permeability(fine_mesh_path, fp, typek=perm_type, export_permfield=export_permfield, update_permfield=update_permfield)
    create_primal_ids(fp, cp, update=bool_export_primal_id)
    export_primal_ids(fine_mesh_path, fp, coarse_mesh_path, export=bool_export_primal_id)
    create_dual_ids(fp, cp, update=bool_export_dual_id, dual_type=my_dual_type)
    export_dual_ids(fine_mesh_path, fp, export=bool_export_dual_id)

    # set_permeability(fp, typek=perm_type, export_permfield=export_permfield)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp, update=update_nodes_weights)
    transm = set_fine_transmissibility_without_bc(
        save_fine_transm_without_bc,
        matrix_path,
        fp,
        fine_transm_without_bc_name
    )
    monotone_transm = set_monotone_transm(
        save_monotone_transm,
        transm,
        w,
        matrix_path,
        monotone_transm_name
    )
    
    resp = set_fine_transmissibility(
        save_fine_transmissibility,
        fp,
        bc,
        matrix_path,
        fine_transmissibility_name
    )

    OR_AMS = get_OR_AMS(fp)
    OP_AMS = get_op(
        save_op,
        fp,
        cp,
        op_name,
        matrix_path,
        op_toget,
        monotone_transm,
        resp,
        level_str
    )
    
    ########################################
    # from packs.multiscale.unstructured.operators.prolongation.amsu import AmsU
    
    # amsu = AmsU()
    # ams_unstructured = Unstructured2DAmsOperator()
    
    # permutation_dict, map_dict = ams_unstructured.get_permutation_matrix_data_2d(
    #     fp[defnames.get_dual_id_name_by_level(1)],
    #     fp['faces']
    # )
    
    # pcorr = amsu.get_correction_function(
    #     OP_AMS,
    #     permutation_dict['G'],
    #     resp['source'],
    #     resp['transmissibility'].tocsc(),
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )
    # coarse_T = OR_AMS*resp['transmissibility']*OP_AMS
    # coarse_q = OR_AMS*(resp['source'] - resp['transmissibility']*pcorr)
    # pc = spsolve(coarse_T.tocsc(), coarse_q)
    # pms_classic = OP_AMS*pc + pcorr
    #####################################
    
    # ##########################################
    # ### iterative ms classic
    # from packs.multiscale.ms_solvers.iterative_solver import iterative_ms_ilu0_bicgstab
    # pit, it = iterative_ms_ilu0_bicgstab(
    #     resp['transmissibility'],
    #     resp['source'],
    #     OP_AMS,
    #     OR_AMS,
    #     epsilon=1e-13
    # )
    # import pdb; pdb.set_trace()
    # ##########################################


    fine_mesh_properties = fp

    if fine_level_setup == 0:
        dual_in_boundary = define_new_fine_levels_v0(fp, bc)
    elif fine_level_setup == 1:
        dual_in_boundary = define_new_fine_levels_v1(fp, bc)
    elif fine_level_setup == 2:
        dual_in_boundary = define_new_fine_levels_v2(fp, bc)
    else:
        raise ValueError
        
    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
    fine_levels = np.full(len(fp['faces']), -1)
    if dual_in_boundary.shape[0] > 0:
        fine_levels[dual_in_boundary] = 0
    
    
    fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
        OR_AMS,
        OP_AMS,
        transm['transmissibility_without_bc'],
        alpha_lim=alpha_lim_finescale
    )

    beta_groups, beta_ind, betas = nu_adm_funcs.get_beta_groups(
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)],
        sp.find(OP_AMS)[0:3],
        fine_mesh_properties['adjacencies'][fine_mesh_properties.internal_edges],
        beta_lim=beta_lim
    )

    # plot_beta_ind(betas, fine_mesh_properties['faces'], 0.3, fine_mesh_path)

    finescale_faces = nu_adm_funcs.get_finescale_vols(
        fine_mesh_properties['faces'][fine_levels==0],
        fine_ids_from_alpha,
        beta_ind,
        beta_groups
    )

    fine_levels[finescale_faces] = 0

    fine_levels[fine_levels==-1] = 1
    finescale_ids = fine_mesh_properties['faces'][fine_levels==0]

    if export_adm_levels_file is True:
        export_adm_levels(fine_mesh_path, fine_levels)
        # flying_mesh_path = _create_flying_mesh(fine_mesh_path)
        test = fine_levels == 0
        if test.sum() > 0:
            print_adm_interfaces_2d(
                fine_mesh_properties,
                fine_mesh_path,
                fine_levels,
                'adm_edges'
            )
        else:
            print("Nao tem volumes na malha fina \n")
    
    LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
        finescale_ids,
        fine_levels,
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)],
        fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    )

    OP_adm, OR_adm = nu_adm_funcs.organize(
        fine_levels,
        sp.find(OP_AMS)[0:3],
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)],
        LEVEL_ID_1,
        ADM_COARSE_ID_LEVEL_1,
        fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    )
    
    ##########################################
    ### iterative ms NU-ADM
    # from packs.multiscale.ms_solvers.iterative_solver import iterative_ms_ilu0_bicgstab
    # pit, it, resid = iterative_ms_ilu0_bicgstab(
    #     resp['transmissibility'],
    #     resp['source'],
    #     OP_AMS,
    #     OR_AMS,
    #     epsilon=1e-13
    # )
    # import pdb; pdb.set_trace()
    ##########################################

    # export_nu_adm_op(fine_mesh_path, OP_adm)

    T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    Q_adm = OR_adm*resp['source']
    P_adm = spsolve(T_adm.tocsc(), Q_adm)
    P_prol = OP_adm*P_adm

    error = np.absolute(pressure - P_prol)
    relative_error = (error/pressure)
    l2_norm_error = np.linalg.norm(error)
    l2_norm_relative_error = np.linalg.norm(relative_error)
    
    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, elements_type='faces')

    mesh_data.create_tag('adm_prol_pressure')
    mesh_data.insert_tag_data('adm_prol_pressure', P_prol, elements_type='faces')
    
    mesh_data.create_tag('absolute_error')
    mesh_data.insert_tag_data('absolute_error', error, elements_type='faces')
    
    mesh_data.create_tag('relative_error')
    mesh_data.insert_tag_data('relative_error', relative_error, elements_type='faces')

    export_adm_name = 'adm_solution_w_' + str(w) + '_' + op_toget + "_NU-ADM"
  
    mesh_data.export_all_elements_type_to_vtk(export_adm_name, element_type='faces')
    
    write_results(
        l2_error=l2_norm_error,
        l2_relative_error=l2_norm_relative_error,
        max_abs_error=error.max(),
        max_abs_relative_error=relative_error.max(),
        dual_type=my_dual_type,
        fine_mesh_name=fp['mesh_name'][0],
        coarse_mesh_name=cp['mesh_name'][0],
        number_fine_vols=fp['faces'].shape[0],
        number_coarse_vols=cp['faces'].shape[0],
        perm_type=perm_type,
        fine_levels_setup=fine_level_setup
    )
    # print('###############################################')
    # print(f' Max abs error {error.max()}')
    # print(f' Max relative error {relative_error.max()}')
    # print(f' L2 error {l2_norm_error}')
    # print(f' L2 relative error {l2_norm_relative_error}')
    # print('###############################################')



