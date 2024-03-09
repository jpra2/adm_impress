from .create_primal_test import get_fine_mesh_path_and_mesh_properties_name_for_test, get_coarse_mesh_path_and_mesh_properties_name_for_test
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio_if_not_exists, _create_flying_mesh
from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d import create_dual
from packs import defnames
from packs.manager import MeshProperty, MeshData

def create_primal_ids(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty):
    if defnames.fine_primal_id not in fine_mesh_properties.keys():
    
        fine_primal_ids = create_coarse_volumes(
            faces_id_level0=fine_mesh_properties['faces'],
            faces_centroids_level0=fine_mesh_properties['faces_centroids'],
            faces_ids_level1=coarse_mesh_properties['faces'],
            nodes_centroids_level1=coarse_mesh_properties['nodes_centroids'],
            nodes_of_faces_level1=coarse_mesh_properties['nodes_of_faces'],
            adjacencies_level0=fine_mesh_properties['adjacencies'],
            faces_of_faces_level0=fine_mesh_properties.faces_of_faces,
            faces_centroids_level1=coarse_mesh_properties['faces_centroids'],
            faces_of_faces_level1=coarse_mesh_properties.faces_of_faces
        )

        fine_mesh_properties.insert_or_update_data({
            defnames.fine_primal_id: fine_primal_ids
        })

        fine_mesh_properties.export_data()



def run():
    fine_mesh_path, fine_mesh_properties_name = get_fine_mesh_path_and_mesh_properties_name_for_test()
    coarse_mesh_path, coarse_mesh_properties_name = get_coarse_mesh_path_and_mesh_properties_name_for_test()

    fine_mesh_properties = create_meshproperties_from_meshio_if_not_exists(fine_mesh_path, fine_mesh_properties_name)
    coarse_mesh_properties = create_meshproperties_from_meshio_if_not_exists(coarse_mesh_path, coarse_mesh_properties_name)

    create_primal_ids(fine_mesh_properties, coarse_mesh_properties)   

    dual_id = create_dual(fine_mesh_properties=fine_mesh_properties, coarse_mesh_properties=coarse_mesh_properties)

    fine_mesh_properties.insert_or_update_data(
        {
            defnames.fine_dual_id: dual_id
        }
    )


    flying_fine_mesh_path = _create_flying_mesh(fine_mesh_path)
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(defnames.fine_dual_id, data_type='int')
    mesh_data.insert_tag_data(defnames.fine_dual_id, dual_id, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements('fine_dual_ids', element_type='faces', elements_array=fine_mesh_properties['faces'])

    # flying_coarse_mesh_path = _create_flying_mesh(coarse_mesh_path)
    # mesh_data = MeshData(mesh_path=flying_coarse_mesh_path)   
    # mesh_data.create_tag('id', data_type='int')
    # mesh_data.insert_tag_data('id', coarse_mesh_properties['faces'], elements_type='faces', elements_array=coarse_mesh_properties['faces'])
    # mesh_data.export_only_the_elements('coarse_ids_test', element_type='faces', elements_array=coarse_mesh_properties['faces'])


    print('fim')



