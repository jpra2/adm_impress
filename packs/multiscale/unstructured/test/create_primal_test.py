from packs import defpaths
import os
# from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio_if_not_exists, _create_flying_mesh
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
from packs import defnames
from packs.manager.mesh_data import MeshData
from packs.utils.multiscale_methods import print_fine_interfaces_coarse_mesh_2d
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh

def get_fine_mesh_path_and_mesh_properties_name_for_test():
    # fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0.msh')
    # fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_2.msh')
    # fine_mesh_path_name_v4 = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_2_v4.msh')
    fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_3.msh')
    fine_mesh_path_name_v4 = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_3_v4.msh')
    # fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_test.msh')
    # fine_mesh_path_name_v4 = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_test_v4.msh')
    # fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0_3.msh')
    fine_mesh_properties_name = 'fine_properties_uns'
    return fine_mesh_path, fine_mesh_properties_name, fine_mesh_path_name_v4

def get_coarse_mesh_path_and_mesh_properties_name_for_test():
    # coarse_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'meshc1_1.msh')
    # coarse_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'meshc1_2.msh')
    coarse_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'meshc3_1.msh')
    coarse_mesh_properties_name = 'coarse_properties_uns'
    return coarse_mesh_path, coarse_mesh_properties_name

def run():
    fine_mesh_path, fine_mesh_properties_name, fine_mesh_path_v4 = get_fine_mesh_path_and_mesh_properties_name_for_test()
    coarse_mesh_path, coarse_mesh_properties_name = get_coarse_mesh_path_and_mesh_properties_name_for_test()

    # fine_mesh_properties = create_meshproperties_from_meshio_if_not_exists(fine_mesh_path, fine_mesh_properties_name)
    # coarse_mesh_properties = create_meshproperties_from_meshio_if_not_exists(coarse_mesh_path, coarse_mesh_properties_name)

    fine_mesh_properties = preprocess_mesh(mesh_name=fine_mesh_path, mesh_properties_name=fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)
    coarse_mesh_properties = preprocess_mesh(mesh_name=coarse_mesh_path, mesh_properties_name=coarse_mesh_properties_name)

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

    key_str = list(fine_primal_ids.keys())[0]
    data = list(fine_primal_ids.values())[0]


    flying_fine_mesh_path = fine_mesh_path
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(key_str, data_type='int')
    mesh_data.insert_tag_data(key_str, data, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements(key_str, element_type='faces', elements_array=fine_mesh_properties['faces'])
    # import pdb; pdb.set_trace()

    print_fine_interfaces_coarse_mesh_2d(
        fine_mesh_properties,
        flying_fine_mesh_path,
        1,
        'edges_selected_2'
    )

    coarse_mesh_data = MeshData(mesh_path=coarse_mesh_path)
    coarse_mesh_data.export_all_elements_type_to_vtk('background_coarse_mesh', element_type='faces')


    
    


