from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility
from packs.multiscale.unstructured.test.test_brazil import define_faces_in_losangle, set_permeability
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation

import os
import numpy as np
from typing import Tuple
from scipy.sparse.linalg import spsolve

def get_properties() -> Tuple[MeshProperty, str]:
    rel_path = os.path.join(
        defpaths.unstructured_coarse_test_mesh_folder,
        'brazil'
    )
    fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
    fine_mesh_properties_name = 'brazilf' 
    fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])

    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    face_p1 = faces[dists <= dists.min()][0]
    dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    face_p0 = faces[dists <= dists.min()][0]
    faces_pressure = np.array([face_p1, face_p0])
    pressure_presc = np.array([1.0, 0.0])

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)
    bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))

    walls_edges = fine_properties['edges'][fine_properties['bool_boundary_edges']]

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
    bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

    return bc

def update_xi_params(xi_params, total_mobility_edges):
    xi_params_new = xi_params.copy()
    xi_params_new[:] = xi_params*total_mobility_edges[:, np.newaxis]
    return xi_params_new
    


def run():

    matrix_path = 'matrices.h5'
    fine_transmissibility_name = 'fine_transmissibility'
    save_fine_transmissibility = True
    update_nodes_weights = False
    type_k = 'barrier'
    export_permfield = False
    update_permfield = False

    relative_perm = BrooksAndCorey()
    biphasic_mobility = BiphasicMobility()
    lsds = LsdsFluxCalculation()

    fp, fine_mesh_path = get_properties()
    bc = set_boundary_conditions(fp)
    define_faces_in_losangle(fp)
    set_permeability(fine_mesh_path, fp, typek=type_k, export_permfield=export_permfield, update_permfield=update_permfield)
    set_weights_nodes(fp, update=update_nodes_weights)
    fp.backup_data('xi_params', 'xi_params_backup')

    fp.insert_or_update_data({
        'saturations': np.repeat(0.2, fp['faces'].shape[0])
    })
    fp.export_data()

    krw_faces, kro_faces = relative_perm.calculate(fp['saturations'])
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    mobw_edges, mobo_edges = direct_edges_mobility(
        mobw_faces,
        mobo_faces,
        fp['areas'],
        fp['faces_of_nodes'],
        fp['nodes_of_edges'],
        bc
    )
    total_mobility_edges = biphasic_mobility.get_total_mobility(mobw_edges, mobo_edges)
    
    fw_edges = biphasic_mobility.get_water_fractionary_flux(mobw_edges, mobo_edges)
    fw_faces = biphasic_mobility.get_water_fractionary_flux(mobw_faces, mobo_faces)

    fp.insert_or_update_data({
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })

    resp = set_fine_transmissibility(
        save_fine_transmissibility,
        fp,
        bc,
        matrix_path,
        fine_transmissibility_name
    )

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

    edges_flux = lsds.get_edges_flux(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    faces_flux = lsds.get_faces_flux(
        edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, 'faces')
    mesh_data.create_tag('faces_flux')
    mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
    mesh_data.export_all_elements_type_to_vtk('pressure_faces', 'faces')




    















    print('fim')
