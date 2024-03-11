from packs.manager.meshmanager import MeshInit
from pymoab import core, types, rng, topo_util
from packs.errors.err import TagNameExistsError, ElementTypeNotInMeshError, DimensionError, MoabTypeNotFoundError
import packs.defpaths as defpaths
import os
import numpy as np
from packs.utils.test_functions import test_mesh_path
from typing import Sequence



class MeshData(MeshInit):
    mesh_elements = ['nodes', 'faces', 'edges', 'volumes']
    
    def __init__(self, dim=2, mesh_path=''):
        self.tags = dict()
        mesh_path_name = test_mesh_path(mesh_path)
        
        self.initialize(mesh_path=mesh_path_name)
        self.init_mesh()
        
        if dim == 2:
            self.init_2d_mesh_entities()
        elif dim == 3:
            self.init_3d_mesh_entities()
        else:
            raise DimensionError('Dimension must be 2 or 3')
    
    def create_tag(self, tag_name, data_size=1, data_type='float'):
        
        self.test_name_in_tags(tag_name)
        
        # if data_density == "dense":
        #     data_density = types.MB_TAG_DENSE
        # elif data_density == "sparse":
        #     data_density = types.MB_TAG_SPARSE
        # elif data_density == "bit":
        #     data_density = types.MB_TAG_BIT
        # else:
        #     print("Please define a valid tag type")
        
        moab_type = self.get_data_type(data_type)
        
        tag = self.mb.tag_get_handle(tag_name, data_size, moab_type, types.MB_TAG_SPARSE, True)

        self.tags.update({
            tag_name: {
                'tag': tag,
                'data_type': data_type,
                'data_size': data_size
            }
        })

    def insert_tag_data(self, tag_name, data, elements_type, elements_array):
        
        all_elements = self.get_all_elements(elements_type)
        
        to_elements = np.array(all_elements).astype(np.uint64)[elements_array]
        
        self.mb.tag_set_data(self.tags[tag_name]['tag'], to_elements, data)
    
    def get_all_elements(self, elements_type):
        
        self.test_elements_type(elements_type)
        
        if elements_type == 'faces':
            all_elements = self.all_faces
        elif elements_type == 'volumes':
            all_elements = self.all_volumes
        elif elements_type == 'edges':
            all_elements = self.all_edges
        elif elements_type == 'nodes':
            all_elements = self.all_nodes
        
        return all_elements
    
    def test_elements_type(self, elements_type):
        
        mesh_elements = self.mesh_elements
        
        if elements_type not in mesh_elements:
            raise ElementTypeNotInMeshError
    
    def test_name_in_tags(self, tag_name):
        
        if not self.verify_tag_in_tags(tag_name):
            pass
        else:
            raise TagNameExistsError(f'The tag {tag_name} already exists')
    
    def verify_tag_in_tags(self, tag_name):
        return tag_name in self.tags

    def export_all_elements_type_to_vtk(self, export_name, element_type):
        
        name = os.path.join(defpaths.results, export_name + '.vtk')
        
        all_elements = self.get_all_elements(element_type)
        meshset = self.mb.create_meshset()
        self.mb.add_entities(meshset, all_elements)
        
        self.mb.write_file(name, [meshset])
    
    def export_only_the_elements(self, export_name, element_type, elements_array):
        
        name = os.path.join(defpaths.results, export_name + '.vtk')
        all_elements = self.get_all_elements(element_type)
        to_elements = np.array(all_elements).astype(np.uint64)[elements_array]
        meshset = self.mb.create_meshset()
        self.mb.add_entities(meshset, to_elements)
        
        self.mb.write_file(name, [meshset])
    
    def get_data_type(self, data_type):
        
        if data_type == 'float':
            moab_type = types.MB_TYPE_DOUBLE
        elif data_type == "int":
            moab_type = types.MB_TYPE_INTEGER
        elif data_type == "bool":
            moab_type = types.MB_TYPE_BIT
        else:
            raise MoabTypeNotFoundError
            
        
        return moab_type
    
    def get_tag_name_for_database(self, tag_name: str, n: int):
        return '_'.join([tag_name, str(n)])
    
    def get_export_name_for_database(self, export_name: str, n: int):
        return '_'.join([export_name, str(n)])

    def insert_array_tag_data(self, tag_name: str, list_data: Sequence[np.ndarray], elements_type: str, list_of_elements_array: Sequence[np.ndarray], data_type='float', data_size=1):

        for i, elements_array in enumerate(list_of_elements_array):
            tag_name_database = self.get_tag_name_for_database(tag_name, i)
            self.create_tag(tag_name=tag_name_database, data_size=data_size, data_type=data_type)
            data = list_data[i]
            assert data.shape[0] == elements_array.shape[0]
            self.insert_tag_data(
                tag_name=tag_name_database,
                data=data,
                elements_type=elements_type,
                elements_array=elements_array
            )

    def export_list_elements_array_data(self, export_name_folder, element_type, list_of_elements_array: Sequence[np.ndarray]):
        folder_to_make = os.path.join(defpaths.results, export_name_folder)
        if not os.path.exists(folder_to_make):
            os.makedirs(folder_to_make)
        
        for i, elements_array in enumerate(list_of_elements_array):
            to_export_name = os.path.join(export_name_folder, '_' + str(i))
            self.export_only_the_elements(
                to_export_name,
                element_type,
                elements_array
            )






            