from packs.multiscale.unstructured.test.create_primal_test import get_fine_mesh_path_and_mesh_properties_name_for_test, get_coarse_mesh_path_and_mesh_properties_name_for_test
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import create_coarse_volumes
# from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio_if_not_exists, _create_flying_mesh
from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d import create_dual
# from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d_v2 import create_dual
from packs import defnames
from packs.manager import MeshProperty, MeshData
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh

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
            level=1,
            edges_ids_level0=fine_mesh_properties['edges'],
            bool_boundary_edges_level0=fine_mesh_properties['bool_boundary_edges'],
            edges_centroids_level0=fine_mesh_properties.edges_centroids
        )

        fine_mesh_properties.insert_or_update_data(
            fine_primal_ids
        )

        fine_mesh_properties.export_data()

def run2():
    fine_mesh_path, fine_mesh_properties_name, fine_mesh_name_v4 = get_fine_mesh_path_and_mesh_properties_name_for_test()
    coarse_mesh_path, coarse_mesh_properties_name = get_coarse_mesh_path_and_mesh_properties_name_for_test()

    # fine_mesh_properties = create_meshproperties_from_meshio_if_not_exists(fine_mesh_path, fine_mesh_properties_name)
    # coarse_mesh_properties = create_meshproperties_from_meshio_if_not_exists(coarse_mesh_path, coarse_mesh_properties_name)

    fine_mesh_properties = preprocess_mesh(mesh_name=fine_mesh_path, mesh_properties_name=fine_mesh_properties_name, mesh_name_v4=fine_mesh_name_v4)
    coarse_mesh_properties = preprocess_mesh(mesh_name=coarse_mesh_path, mesh_properties_name=coarse_mesh_properties_name)

    create_primal_ids(fine_mesh_properties, coarse_mesh_properties)   

    dual_data = create_dual(fine_mesh_properties=fine_mesh_properties, coarse_mesh_properties=coarse_mesh_properties, level=1)

    fine_mesh_properties.insert_or_update_data(
        dual_data
    )

    fine_mesh_properties.export_data()

    key_str = defnames.get_dual_id_name_by_level(level=1)
    data = dual_data.get(key_str)

    flying_fine_mesh_path = fine_mesh_path
    mesh_data = MeshData(mesh_path=flying_fine_mesh_path)   
    mesh_data.create_tag(key_str, data_type='int')
    mesh_data.insert_tag_data(key_str, data, elements_type='faces', elements_array=fine_mesh_properties['faces'])
    mesh_data.export_only_the_elements(key_str, element_type='faces', elements_array=fine_mesh_properties['faces'])

    dual_volumes_name = defnames.get_dual_volumes_name_by_level(1)
    dual_volumes = dual_data[dual_volumes_name]

    mesh_data.export_list_elements_array_data(dual_volumes_name, 'faces', dual_volumes)

    interaction_regions_name = defnames.get_dual_interation_region_name_by_level(1)
    regions = dual_data[interaction_regions_name]

    mesh_data.export_list_elements_array_data(interaction_regions_name, 'faces', regions)

    dual_boundarys_name = 'dual_boundary'
    boundarys = dual_data[defnames.boundary_dual_interaction + defnames.level_str(1)]

    mesh_data.export_list_elements_array_data(dual_boundarys_name, 'faces', boundarys)

    # selected_edges_name = defnames.edges_selected + '_level' + str(1) 
    # edges_data = dual_data[selected_edges_name]
    # mesh_data.export_list_elements_array_data(selected_edges_name, 'faces', edges_data)




    # flying_coarse_mesh_path = _create_flying_mesh(coarse_mesh_path)
    # mesh_data = MeshData(mesh_path=flying_coarse_mesh_path)   
    # mesh_data.create_tag('id', data_type='int')
    # mesh_data.insert_tag_data('id', coarse_mesh_properties['faces'], elements_type='faces', elements_array=coarse_mesh_properties['faces'])
    # mesh_data.export_only_the_elements('coarse_ids_test', element_type='faces', elements_array=coarse_mesh_properties['faces'])


    print('fim')



