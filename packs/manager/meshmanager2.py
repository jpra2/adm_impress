import numpy as np
import pandas as pd
import os
import copy

from packs.manager.arraydatamanager import ArrayDataManager, test_str_instance, test_array_instance
from packs.errors import err as errors
from packs import defpaths

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
            return os.path.join(defpaths.data_mesh, 'mesh_property_' + self.mesh_name[0] + '.npz')
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
            if self.verify_name_in_data_names(name):
                del self.__dict__[name]

    def backup_data(self, from_name: str, to_name: str):
        
        self.verify_name_in_data_names_or_raise_error(from_name)
        data = self[from_name].copy()
        self.insert_or_update_data({
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

    @property
    def edges_dim(self):
        try:
            return self['edges_dim']
        except KeyError:
            resp = np.linalg.norm(
                self.nodes_centroids[self.nodes_of_edges[self.edges, 0]] - self.nodes_centroids[self.nodes_of_edges[self.edges, 1]],
                axis=1
            )
            self.insert_data({'edges_dim': resp})
            self.export_data()
            return resp
    
    @property
    def edges_centroids(self):
        try:
            return self['edges_centroids']
        except KeyError:
            resp = (self.nodes_centroids[self.nodes_of_edges[:, 1]] + self.nodes_centroids[self.nodes_of_edges[:, 0]])/2
            self.insert_data({'edges_centroids': resp})
            self.export_data()
            return resp
        # resp = np.mean(
        #     self.nodes_centroids[self.nodes_of_edges],
        #     axis=1
        # )
    
    @property
    def boundary_edges(self):
        return self.edges[self.bool_boundary_edges]
    
    @property
    def internal_edges(self):
        return np.setdiff1d(self.edges, self.boundary_edges, assume_unique=True)
    
    @property
    def boundary_nodes(self):
        return self.nodes[self.bool_boundary_nodes]
    
    @property
    def internal_nodes(self):
        return np.setdiff1d(self.nodes, self.boundary_nodes, assume_unique=True)


    @property
    def faces_of_faces(self):
        try:
            return self['faces_of_faces']
        except KeyError:
            faces_of_faces = []
            for face in self.faces:
                test1 = self.adjacencies[:, 0] == face
                test2 = self.adjacencies[:, 1] == face
                test3 = test1 | test2
                adjs = self.adjacencies[test3]
                adjs = adjs[adjs!=face]
                adjs = adjs[adjs!=-1]
                faces_of_faces.append(adjs)
            
            faces_of_faces = np.array(faces_of_faces, dtype='O')
            self.insert_data({'faces_of_faces': faces_of_faces})
            self.export_data()
            return faces_of_faces
        
    @property
    def faces_of_faces_by_nodes(self):
        try:
            return self['faces_of_faces_by_nodes']
        except KeyError:
            faces_of_faces_by_nodes = []
            for face in self.faces:
                nodes_face = self['nodes_of_faces'][face]
                faces_nodes_face = np.concatenate(
                    self['faces_of_nodes'][nodes_face]
                )
                faces_nodes_face = np.setdiff1d(faces_nodes_face, [face])
                faces_of_faces_by_nodes.append(faces_nodes_face)
            
            faces_of_faces_by_nodes = np.array(faces_of_faces_by_nodes, dtype='O')
            self.insert_data({'faces_of_faces_by_nodes': faces_of_faces_by_nodes})
            self.export_data()
            return faces_of_faces_by_nodes
    
    def get_nodes_org_from_faces_of_nodes_object(self):
        n_faces_of_nodes = self['n_faces_of_nodes']
        faces_of_nodes = self['faces_of_nodes']

        all_nodes_org = []
        faces_of_nodes_org = []
        nf_nodes = np.unique(n_faces_of_nodes)

        for i in nf_nodes:
            test = n_faces_of_nodes == i
            v4 = self.nodes[test]
            all_nodes_org.append(v4)
        
            ft = faces_of_nodes[v4].copy()
            ft2 = np.concatenate(ft)
            ft3 = ft2.reshape((ft.shape[0], i))
            faces_of_nodes_org.append(ft3)
        
        all_nodes_org = np.array(all_nodes_org, dtype='O')
        faces_of_nodes_org = np.array(faces_of_nodes_org, dtype='O')

        resp = {
            'nodes_org': all_nodes_org,
            'faces_of_nodes_org': faces_of_nodes_org,
            'n_nodes_org': nf_nodes
        }

        return resp

    def get_internal_nodes_org_from_faces_of_nodes_object(self):
        if self.verify_name_in_data_names('internal_nodes_org'):
            return {
                'internal_nodes_org': self['internal_nodes_org'],
                'internal_faces_of_nodes_org': self['internal_faces_of_nodes_org'],
                'internal_n_nodes': self['internal_n_nodes']
            }   
        n_faces_of_nodes = self['n_faces_of_nodes']
        internal_nodes = self.internal_nodes
        faces_of_nodes = self['faces_of_nodes']

        all_nodes_org = []
        faces_of_nodes_org = []
        n_nodes = np.arange(2, n_faces_of_nodes.max()+1)
        new_n_nodes = []

        for i in n_nodes:
            test = n_faces_of_nodes == i
            v4 = self.nodes[test]
            v4 = np.intersect1d(v4, internal_nodes)
            if v4.shape[0] == 0:
                continue
            all_nodes_org.append(v4)
        
            ft = faces_of_nodes[v4].copy()
            ft2 = np.concatenate(ft)
            ft3 = ft2.reshape((ft.shape[0], i))
            faces_of_nodes_org.append(ft3)
            new_n_nodes.append(i)
        
        all_nodes_org = np.array(all_nodes_org, dtype='O')
        faces_of_nodes_org = np.array(faces_of_nodes_org, dtype='O')
        new_n_nodes = np.array(new_n_nodes)

        resp = {
            'internal_nodes_org': all_nodes_org,
            'internal_faces_of_nodes_org': faces_of_nodes_org,
            'internal_n_nodes': new_n_nodes
        }
        
        return resp

    @property
    def dist_centroids(self):
        try:
            return self['dist_centroids']
        except KeyError:
            all_dists = np.zeros(len(self['edges']))
            faces_centroids = self['faces_centroids']
            adjacencies = self['adjacencies']
            bool_internal_edges = ~self['bool_boundary_edges']
            edges_centroids = self.edges_centroids

            all_dists[self.internal_edges] = np.linalg.norm(
                faces_centroids[adjacencies[bool_internal_edges, 0]] - faces_centroids[adjacencies[bool_internal_edges, 1]],
                axis=1
            )

            all_dists[self.boundary_edges] = np.linalg.norm(
                faces_centroids[adjacencies[self.boundary_edges, 0]] - edges_centroids[self.boundary_edges],
                axis=1
            )

            self.insert_data({'dist_centroids': all_dists})
            self.export_data()
            return all_dists

    @property
    def faces(self):
        return self['faces']
    
    @property
    def nodes(self):
        return self['nodes']
    
    @property
    def edges(self):
        return self['edges']
    
    @property
    def bedges_without_bedges_to_remove(self):
        edges_to_remove = self.remove_bedges
        
        if edges_to_remove.shape[0] == 0:
            return self.boundary_edges
        else:
            my_edges = np.setdiff1d(self.boundary_edges, edges_to_remove)
            return my_edges
        
    @property
    def remove_bedges(self):
        try:
            edges_to_remove = self['remove_bedges']
        except KeyError:
            edges_to_remove = np.array([])
        
        return edges_to_remove