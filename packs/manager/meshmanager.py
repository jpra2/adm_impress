# from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
import numpy as np
import pandas as pd
import copy
import packs.defpaths as defpaths
import os
from packs.manager.arraydatamanager import ArrayDataManager, test_str_instance, test_array_instance
from packs.errors import err as errors
from packs.utils import calculate_face_properties
import copy


class MeshInit:
    
    def initialize(self, mesh_path='', mesh_name=''):

        self.mesh_path = mesh_path
        self.all_volumes = None
        self.all_faces = None
        self.all_nodes = None
        self.all_edges = None
        self.mb: core.Core = None
        self.mtu: topo_util.MeshTopoUtil = None
        self.root_set = None
        self.data = dict()
        self.mesh_name = mesh_name

    def _init_mesh(self):
        mb = core.Core()
        scd = ScdInterface(mb)
        mtu = topo_util.MeshTopoUtil(mb)
        try:
            mb.load_file(self.mesh_path)
        except OSError:
            raise OSError(f'Erro na leitura da malha: {self.mesh_path}')
        root_set = mb.get_root_set()
        return mb, mtu, root_set

    def init_mesh(self):
        self.mb, self.mtu, self.root_set = self._init_mesh()

    def init_3d_mesh_entities(self):
        self.all_volumes = self.mb.get_entities_by_dimension(0, 3)
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        boundary_faces = self.mb.get_entities_by_dimension(0, 2)
        # edges = self.mb.get_entities_by_dimension(0, 1)
        
        # all_entities = self.mb.get_entities_by_handle(self.root_set)
        
        # all_entities = np.setdiff1d(all_entities, self.all_volumes)
        # all_entities = np.setdiff1d(all_entities, self.all_nodes)
        # all_entities = np.setdiff1d(all_entities, faces)
        # all_entities = np.setdiff1d(all_entities, edges)
        
        # for ent in all_entities:
        #     ents = self.mb.get_entities_by_handle(ent)
        #     print(ents)
        
        
        # type_moab = self.mb.type_from_handle(self.root_set)
        all_moab_types = dir(types)
        
        self.dict_moab_types = dict()
        
        for tt in all_moab_types[0:-20]:
            
            # exec('print(types.' + tt + ')')
            # exec('respp[tt] = types.' + tt)
            exec('self.dict_moab_types[types.' + tt +'] = tt')
        
        
        self.mtu.construct_aentities(self.all_nodes)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)
        # other_faces = []
        # for face in self.all_faces:
        #     tags = self.mb.tag_get_tags_on_entity(face)
        #     if len(tags) != 1:
        #         other_faces.append(face)
        
        # boundary_faces = np.setdiff1d(self.all_faces, other_faces)
        
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)

    def init_2d_mesh_entities(self):
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)
        
        # type_moab = self.mb.type_from_handle(self.root_set)
        all_moab_types = dir(types)
        
        self.dict_moab_types = dict()
        
        for tt in all_moab_types[0:-20]:
            
            # exec('print(types.' + tt + ')')
            # exec('respp[tt] = types.' + tt)
            exec('self.dict_moab_types[types.' + tt +'] = tt')
        
        boundary_edges = self.mb.get_entities_by_dimension(0, 1)
        self.mtu.construct_aentities(self.all_nodes)
        
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)
    


class MeshProperty:
    
    def insert_mesh_name(self, name=''):
        self.__dict__['mesh_name'] = np.array([name])
        
    
    def insert_data(self, data: dict):
        """data is a dictionary with str keys and np.ndarray values

        Args:
            data (_type_): dict
        """
        names = list(data.keys())
        values = list(data.values())
        
        a = [test_array_instance(_) for _ in values]
        a = [test_str_instance(_) for _ in names]
        
        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(list(self.__dict__.keys()))
        names_data_self = names_data_self[names_data_self != 'mesh_name']
        test = names_series.isin(names_data_self)
        if test.any().values[0]:
            names_in = names_series[test.values].values.flatten()
            raise errors.NameExistsError(f'The names: - {names_in} - exists in mesh properties')
        
        
        self.__dict__.update(data)
        # self._data.update(data)
        # self.__dict__['mesh_name'] = data['mesh_name'][0]
    
    def __setattr__(self, name, value):
        raise Exception("It is read only!")      
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    @property
    def class_path(self):
        try:
            return os.path.join(defpaths.flying, 'mesh_property_' + self.mesh_name[0] + '.npz')
        except:
            import pdb; pdb.set_trace()
    
    def export_data(self):
        manager = ArrayDataManager(self.class_path)
        manager.insert_data(self.__dict__)
        manager.export()

    def load_data(self):
        self.verify_if_exists()
        
        manager = ArrayDataManager(self.class_path)
        self.insert_data(manager.get_data_from_load())

    def get_all_data(self):
        return self.__dict__

    def rename_data(self, datas_to_rename: dict):
        """ Update the data name

        @param datas_to_rename: dict where key = old data name, value = new data name
        """

        new_data = dict()

        for name in list(datas_to_rename.keys()):
            self.verify_name_in_data_names_or_raise_error(name)
            self.verify_name_not_in_data_names_or_raise_error(datas_to_rename[name])
            data = copy.deepcopy(self[name])
            new_name = datas_to_rename[name]
            del self.__dict__[name]
            new_data.update({new_name: data})

        self.insert_data(new_data)
    
    def update_data(self, datas_to_update: dict):
        
        new_data = dict()
        for name in datas_to_update:
            self.verify_name_in_data_names_or_raise_error(name)
            new_data.update({
                name: datas_to_update[name]   
            })
            del self.__dict__[name]
        
        self.insert_data(new_data)
    
    def insert_or_update_data(self, datas: dict):
        
        data_names = list(datas.keys())
        names_in = self.verify_names_in_data_names(data_names)
        names_out = self.verify_names_out_data_names(data_names)
        
        self.insert_data({name: datas[name] for name in names_out})
        self.update_data({name: datas[name] for name in names_in})
        
    def remove_data(self, data_name: list):
        for name in data_name:
            del self.__dict__[name]

    def backup_data(self, from_name: str, to_name: str):
        
        self.verify_name_in_data_names_or_raise_error(from_name)
        self.verify_name_not_in_data_names_or_raise_error(to_name)
        data = self[from_name].copy()
        self.insert_data({
            to_name: data
        })
    
    def backup_datas(self, backup_datas_name: dict):

        """
        backup_datas_name = {from_name1: to_name1, from_name2: to_name2 ...}
        """

        for name in backup_datas_name:
            self.backup_data(name, backup_datas_name[name])

    def exists(self):
        return os.path.exists(self.class_path)
    
    def verify_if_exists(self):
        if self.exists():
            pass
        else:
            raise FileExistsError
    
    def keys(self):
        return self.__dict__.keys()
    
    @property        
    def data_names(self):
        return list(self.keys())
    
    def verify_names_in_data_names(self, names: list):
        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(self.data_names)
        test = names_series.isin(names_data_self)
        test = test.values
        names_in = names_series[test].values.flatten()
        return names_in
    
    def verify_names_out_data_names(self, names:list):
        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(self.data_names)
        test = names_series.isin(names_data_self)
        test = ~test.values
        names_out = names_series[test].values.flatten()
        return names_out

    def verify_name_in_data_names(self, name: str):
        return name in self.data_names   

    def verify_name_in_data_names_or_raise_error(self, name: str):
        if self.verify_name_in_data_names(name):
            pass
        else:
            raise errors.NameExistsError(f'The name: - {name} - does not exists in mesh properties')
    
    def verify_name_not_in_data_names_or_raise_error(self, name: str):
        if self.verify_name_in_data_names(name):
            raise errors.NameExistsError(f'The name: - {name} - exists in mesh properties')

        
    

class CreateMeshProperties(MeshInit):

    '''
        Create mesh properties using pymoab
    '''

    def _init_3d_properties(self, faces, volumes, nodes):
        n_faces = len(faces)
        volumes_adj_by_faces = np.repeat(-1, n_faces*2).reshape((n_faces, 2)).astype(np.uint64)
        volumes_series = pd.DataFrame({
            'vol_ids': volumes
        }, index = np.array(self.all_volumes))

        nodes_series = pd.DataFrame({
            'nodes_ids': nodes
        }, index=self.all_nodes)
        
        faces_series = pd.DataFrame({
            'faces_ids': faces
        }, index=self.all_faces)

        # nodes_of_faces = np.repeat(-1, n_faces*4).reshape((n_faces, 4)).astype(np.uint64)
        nodes_of_faces = []
        faces_of_volumes = []
        nodes_of_volumes = []
        volumes_adj_by_nodes = []

        for i, face in enumerate(self.all_faces):
            volumes_adj_by_faces[i][:] = self.mtu.get_bridge_adjacencies(face, 2, 3)
            nodes_of_faces_elems = self.mtu.get_bridge_adjacencies(face, 2, 0) 
            nodes_of_faces.append(nodes_series.loc[nodes_of_faces_elems].to_numpy().flatten())
            
        for vol in self.all_volumes:
            faces_of_volumes_elements = self.mtu.get_bridge_adjacencies(vol, 3, 2)
            nodes_of_volumes_elements = self.mtu.get_bridge_adjacencies(vol, 3, 0)
            vols_by_nodes_elements = self.mtu.get_bridge_adjacencies(vol, 0, 3)
            faces_of_volumes_loc = faces_series.loc[faces_of_volumes_elements].to_numpy().flatten()
            nodes_of_volumes_loc = nodes_series.loc[nodes_of_volumes_elements].to_numpy().flatten()
            vols_by_nodes_loc = volumes_series.loc[vols_by_nodes_elements].to_numpy().flatten()
            faces_of_volumes.append(faces_of_volumes_loc)
            nodes_of_volumes.append(nodes_of_volumes_loc)
            volumes_adj_by_nodes.append(vols_by_nodes_loc)
        

        test = volumes_adj_by_faces[:, 0] == volumes_adj_by_faces[:, 1]
        volumes_adj_by_faces[:, 0] = volumes_series.loc[volumes_adj_by_faces[:, 0]].to_numpy().flatten()
        volumes_adj_by_faces[:, 1] = volumes_series.loc[volumes_adj_by_faces[:, 1]].to_numpy().flatten()
        volumes_adj_by_faces = volumes_adj_by_faces.astype(np.int64)
        volumes_adj_by_faces[test, 1] = -1
        
        nodes_of_faces = np.array(nodes_of_faces)
        faces_of_volumes = np.array(faces_of_volumes)
        nodes_of_volumes = np.array(nodes_of_volumes)
        volumes_adj_by_nodes = np.array(volumes_adj_by_nodes)
        
        bool_internal_faces = test

        return bool_internal_faces, volumes_adj_by_faces, nodes_of_faces, faces_of_volumes, nodes_of_volumes, volumes_adj_by_nodes
    
    def _init_2d_properties(self, faces, edges, nodes, nodes_centroids):
        n_edges = len(edges)
        faces_adj_by_edges = np.repeat(-1, n_edges*2).reshape((n_edges, 2)).astype(np.uint64)
        nodes_of_edges = faces_adj_by_edges.copy()
        faces_centroids = np.zeros((len(faces), nodes_centroids.shape[1]))

        nodes_series = pd.DataFrame({
            'nodes_ids': nodes
        }, index=self.all_nodes)
        
        faces_series = pd.DataFrame({
            'faces_ids': faces
        }, index=self.all_faces)
        
        edges_series = pd.DataFrame({
            'edges_ids': edges
        }, index=self.all_edges)
        
        # nodes_of_faces = np.repeat(-1, n_faces*4).reshape((n_faces, 4)).astype(np.uint64)
        nodes_of_faces = []
        edges_of_faces = []
        faces_adj_by_nodes = []

        for i, edge in enumerate(self.all_edges):
            faces_adj_by_edges[i][:] = self.mtu.get_bridge_adjacencies(edge, 1, 2)
            nodes_of_edge_elems = self.mtu.get_bridge_adjacencies(edge, 1, 0) 
            nodes_of_edges[i][:]= nodes_series.loc[nodes_of_edge_elems].to_numpy().flatten()
            
        for i, face in enumerate(self.all_faces):
            edges_of_face_elements = self.mtu.get_bridge_adjacencies(face, 2, 1)
            nodes_of_face_elements = self.mtu.get_bridge_adjacencies(face, 2, 0)
            faces_by_nodes_elements = self.mtu.get_bridge_adjacencies(face, 0, 2)
            
            edges_of_faces_loc = edges_series.loc[edges_of_face_elements].to_numpy().flatten()
            
            nodes_of_face_loc = nodes_series.loc[nodes_of_face_elements].to_numpy().flatten()
            
            faces_by_nodes_loc = faces_series.loc[faces_by_nodes_elements].to_numpy().flatten()
            
            edges_of_faces.append(edges_of_faces_loc)
            nodes_of_faces.append(nodes_of_face_loc)
            faces_adj_by_nodes.append(faces_by_nodes_loc)
            faces_centroids[i][:] = np.mean(nodes_centroids[nodes_of_face_loc], axis=0)
        

        test = faces_adj_by_edges[:, 0] == faces_adj_by_edges[:, 1]
        faces_adj_by_edges[:, 0] = faces_series.loc[faces_adj_by_edges[:, 0]].to_numpy().flatten()
        faces_adj_by_edges[:, 1] = faces_series.loc[faces_adj_by_edges[:, 1]].to_numpy().flatten()
        faces_adj_by_edges = faces_adj_by_edges.astype(np.int64)
        faces_adj_by_edges[test, 1] = -1
        
        nodes_of_faces = np.array(nodes_of_faces)
        edges_of_faces = np.array(edges_of_faces)
        faces_adj_by_nodes = np.array(faces_adj_by_nodes, dtype='O')
        bool_boundary_edges = test        
        
        return bool_boundary_edges, faces_adj_by_edges, nodes_of_faces, edges_of_faces, faces_adj_by_nodes, nodes_of_edges, faces_centroids
        
    def get_nodes_and_edges_and_faces_adjacencies_by_nodes(self, nodes, edges, nodes_of_edges, faces):
        
        nodes_series = pd.DataFrame({
            'nodes_ids': nodes
        }, index=self.all_nodes)
        
        edges_series = pd.DataFrame({
            'edges_ids': edges
        }, index=self.all_edges)
        
        faces_series = pd.DataFrame({
            'faces_ids': faces
        }, index=self.all_faces)
        
        nodes_adj_by_nodes = []
        edges_adj_by_nodes = []
        faces_adj_by_nodes = []
        
        for node_elem in self.all_nodes:
            nodes_adj_elem = self.mtu.get_bridge_adjacencies(node_elem, 1, 0)
            edges_adj_elem = self.mtu.get_bridge_adjacencies(node_elem, 0, 1)
            faces_adj_elem = self.mtu.get_bridge_adjacencies(node_elem, 0, 2)
            
            nodes_adj = nodes_series.loc[nodes_adj_elem].to_numpy().flatten()
            edges_adj = edges_series.loc[edges_adj_elem].to_numpy().flatten()
            faces_adj = faces_series.loc[faces_adj_elem].to_numpy().flatten()
            
            node = nodes_series.loc[node_elem].values[0]
            nodes_of_edges_adj = nodes_of_edges[edges_adj]
            
            test = nodes_of_edges_adj[nodes_of_edges_adj != node]
            
            order = pd.DataFrame({
                'index_edges': np.arange(len(edges_adj))
            }, index=test)
            
            order_correct = order.loc[nodes_adj].values.flatten()
            
            nodes_adj_by_nodes.append(nodes_adj)
            edges_adj_by_nodes.append(edges_adj[order_correct])
            faces_adj_by_nodes.append(faces_adj)
            
        
        nodes_adj_by_nodes = np.array(nodes_adj_by_nodes, dtype='O')
        edges_adj_by_nodes = np.array(edges_adj_by_nodes, dtype='O')
        faces_adj_by_nodes = np.array(faces_adj_by_nodes, dtype='O')
        
        return nodes_adj_by_nodes, edges_adj_by_nodes, faces_adj_by_nodes
        
    def create_3d_initial_array_properties(self):

        volumes = np.arange(len(self.all_volumes), dtype=int)
        faces = np.arange(len(self.all_faces), dtype=int)
        edges = np.arange(len(self.all_edges), dtype=int)
        nodes = np.arange(len(self.all_nodes), dtype=int)
        bool_internal_faces, volumes_adj_by_faces, nodes_of_faces, faces_of_volumes, nodes_of_volumes, volumes_adj_by_nodes =  self._init_3d_properties(faces, volumes, nodes)

        nodes_centroids = np.array([self.mb.get_coords(node) for node in self.all_nodes])

        # av2 = np.array(self.all_volumes, np.uint64)
        #
        # vc = self.mtu.get_average_position(av2[0])

        self.data['volumes'] = volumes
        self.data['faces'] = faces
        self.data['edges'] = edges
        self.data['nodes'] = nodes
        self.data['bool_internal_faces'] = bool_internal_faces
        self.data['volumes_adj_by_faces'] = volumes_adj_by_faces
        self.data['volumes_adj_by_nodes'] = volumes_adj_by_nodes
        self.data['nodes_of_faces'] = nodes_of_faces
        self.data['nodes_centroids'] = nodes_centroids
        self.data['faces_of_volumes'] = faces_of_volumes
        self.data['nodes_of_volumes'] = nodes_of_volumes
        self.data['mesh_name'] = np.array([self.mesh_name])
    
    def create_2d_initial_array_properties(self):
        
        faces = np.arange(len(self.all_faces), dtype=int)
        edges = np.arange(len(self.all_edges), dtype=int)
        nodes = np.arange(len(self.all_nodes), dtype=int)
        nodes_centroids = np.array([self.mb.get_coords(node) for node in self.all_nodes])
        
        bool_boundary_edges, faces_adj_by_edges, nodes_of_faces, edges_of_faces, faces_adj_by_nodes, nodes_of_edges, faces_centroids =  self._init_2d_properties(faces, edges, nodes, nodes_centroids)
        
        nodes_adj_by_nodes, edges_adj_by_nodes, faces_adj_by_nodes = self.get_nodes_and_edges_and_faces_adjacencies_by_nodes(nodes, edges, nodes_of_edges, faces)
        
        calculate_face_properties.ordenate_edges_and_nodes_of_nodes_xy_plane(
            nodes,
            edges, 
            nodes_adj_by_nodes,
            edges_adj_by_nodes,
            nodes_centroids
        )
        
        calculate_face_properties.ordenate_faces_of_nodes_xy_plane(faces_centroids, faces_adj_by_nodes, nodes_centroids)
        
        bool_boundary_nodes = calculate_face_properties.define_bool_boundary_nodes(bool_boundary_edges, nodes_of_edges, nodes)
        
        # av2 = np.array(self.all_volumes, np.uint64)
        #
        # vc = self.mtu.get_average_position(av2[0])

        self.data['faces'] = faces
        self.data['edges'] = edges
        self.data['nodes'] = nodes
        self.data['bool_boundary_edges'] = bool_boundary_edges
        self.data['faces_adj_by_edges'] = faces_adj_by_edges
        self.data['faces_adj_by_nodes'] = faces_adj_by_nodes
        self.data['nodes_of_faces'] = nodes_of_faces
        self.data['nodes_centroids'] = nodes_centroids
        self.data['edges_of_faces'] = edges_of_faces
        self.data['mesh_name'] = np.array([self.mesh_name])
        self.data['nodes_of_edges'] = nodes_of_edges
        self.data['faces_centroids'] = faces_centroids
        self.data['nodes_adj_by_nodes'] = nodes_adj_by_nodes
        self.data['edges_adj_by_nodes'] = edges_adj_by_nodes
        self.data['bool_boundary_nodes'] = bool_boundary_nodes
        
    def export_data_msh(self):
        # import pdb; pdb.set_trace()
        tags = self.mb.tag_get_tags_on_entity(self.root_set)
        import pdb; pdb.set_trace()
        
        self.mb.write_file(os.path.join(defpaths.mesh, self.mesh_name) + '.msh')
          
    def create_3d_mesh_data(self):
        self.init_mesh()
        self.init_3d_mesh_entities()
        self.create_3d_initial_array_properties()
        mesh_property = MeshProperty()
        mesh_property.insert_mesh_name(self.data['mesh_name'][0])
        mesh_property.insert_data(copy.deepcopy(self.data))
        mesh_property.export_data()
        return mesh_property
        
    def create_2d_mesh_data(self):
        self.init_mesh()
        self.init_2d_mesh_entities()
        self.create_2d_initial_array_properties()
        mesh_property = MeshProperty()
        mesh_property.insert_mesh_name(self.data['mesh_name'][0])
        mesh_property.insert_data(copy.deepcopy(self.data))
        # mesh_property.export_data()
        return mesh_property
        


def create_initial_mesh_properties(mesh_path, mesh_name):
    mesh_create = CreateMeshProperties()
    mesh_create.initialize(mesh_path=mesh_path, mesh_name=mesh_name)
    mesh_properties: MeshProperty = mesh_create.create_2d_mesh_data()
    return mesh_properties

def load_mesh_properties(mesh_name):
    mesh_properties = MeshProperty()
    mesh_properties.insert_mesh_name(mesh_name)
    mesh_properties.load_data()
    return mesh_properties


# if __name__ == '__main__':
#     mesh_path = 'mesh/80x80x1_ufce.h5m'
#     mesh_name = '80x80_ufce'
#     mesh_create = CreateMeshProperties()
#     mesh_create.initialize(mesh_path=mesh_path, mesh_name=mesh_name)
#     mesh_properties = mesh_create.create_mesh_data()
    
#     import pdb; pdb.set_trace()