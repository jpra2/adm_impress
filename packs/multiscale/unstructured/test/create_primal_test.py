from packs import defpaths
import os
from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio_if_not_exists, _create_flying_mesh
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
from packs import defnames
from packs.manager.mesh_data import MeshData

def get_fine_mesh_path_and_mesh_properties_name_for_test():
    fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0.msh')
    fine_mesh_properties_name = 'fine_properties_uns'
    return fine_mesh_path, fine_mesh_properties_name

def get_coarse_mesh_path_and_mesh_properties_name_for_test():
    coarse_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'meshc1_1.msh')
    coarse_mesh_properties_name = 'coarse_properties_uns'
    return coarse_mesh_path, coarse_mesh_properties_name

def run():
    fine_mesh_path, fine_mesh_properties_name = get_fine_mesh_path_and_mesh_properties_name_for_test()
    coarse_mesh_path, coarse_mesh_properties_name = get_coarse_mesh_path_and_mesh_properties_name_for_test()

    fine_mesh_properties = create_meshproperties_from_meshio_if_not_exists(fine_mesh_path, fine_mesh_properties_name)
    coarse_mesh_properties = create_meshproperties_from_meshio_if_not_exists(coarse_mesh_path, coarse_mesh_properties_name)

    fine_primal_ids = create_coarse_volumes(
        faces_id_level0=fine_mesh_properties['faces'],
        faces_centroids_level0=fine_mesh_properties['faces_centroids'],
        faces_ids_level1=coarse_mesh_properties['faces'],
        nodes_centroids_level1=coarse_mesh_properties['nodes_centroids'],
        nodes_of_faces_level1=coarse_mesh_properties['nodes_of_faces'],
        faces_of_nodes_level0=fine_mesh_properties['faces_of_nodes'],
        adjacencies_level0=fine_mesh_properties['adjacencies'],
        faces_of_faces_level0=fine_mesh_properties.faces_of_faces,
        faces_centroids_level1=coarse_mesh_properties['faces_centroids'],
        faces_of_faces_level1=coarse_mesh_properties.faces_of_faces
    )
    
    fine_mesh_properties.insert_or_update_data(
        {
            defnames.fine_primal_id: fine_primal_ids
        }
    )


    flying_fine_mesh_path = _create_flying_mesh(fine_mesh_path)
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(defnames.fine_primal_id, data_type='int')
    mesh_data.insert_tag_data(defnames.fine_primal_id, fine_primal_ids, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements('fine_prinal_ids', element_type='faces', elements_array=fine_mesh_properties['faces'])
    print('end')
    # import pdb; pdb.set_trace()
    
    


