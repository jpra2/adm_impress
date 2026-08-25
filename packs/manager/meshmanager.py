# from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
# from pymoab import core, types, rng, topo_util
# from pymoab import skinner as sk
# from pymoab.scd import ScdInterface
# from pymoab.hcoord import HomCoord
import numpy as np
import pandas as pd
import copy
import packs.defpaths as defpaths
import os
from packs.manager.arraydatamanager import ArrayDataManager, test_str_instance, test_array_instance
from packs.errors import err as errors
from packs.utils import calculate_face_properties
import copy
from packs.manager.meshmanager2 import MeshProperty

class MeshInit:
    
    def initialize(self, mesh_path='', mesh_name=''):
        
        from pymoab import core, topo_util

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
        from pymoab import core, topo_util
        from pymoab.scd import ScdInterface
        
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
        from pymoab import types
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
        from pymoab import types
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)


        # bfaces = self.mb.get_entities_by_type_and_tag(self.root_set, types.MBTRI)
        
        # type_moab = self.mb.type_from_handle(self.root_set)
        all_moab_types = dir(types)
        
        self.dict_moab_types = dict()
        
        for tt in all_moab_types[0:-20]:
            
            # exec('print(types.' + tt + ')')
            # exec('respp[tt] = types.' + tt)
            exec('self.dict_moab_types[types.' + tt +'] = tt')
        
        # import pdb; pdb.set_trace()
        
        boundary_edges = self.mb.get_entities_by_dimension(0, 1)
        self.mtu.construct_aentities(self.all_nodes)
        
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)

        # select_nodes = self.test_nodes_out_2d()
        # self.all_nodes = np.array(self.all_nodes[select_nodes])


    def test_nodes_out_2d(self):
        ntotalnodes = len(self.all_nodes)
        new_nodes_bool = np.full(ntotalnodes, True, bool)

        for i, node in enumerate(self.all_nodes):
            faces_node = self.mtu.get_bridge_adjacencies(node, 0, 2)
            if len(faces_node) == 0:
                new_nodes_bool[i] = False
        
        return new_nodes_bool
     

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
        
        nodes_adj_by_nodes, edges_adj_by_nodes = calculate_face_properties.ordenate_edges_and_nodes_of_nodes_xy_plane(
            nodes,
            edges, 
            nodes_adj_by_nodes,
            edges_adj_by_nodes,
            nodes_centroids
        )
        
        faces_adj_by_nodes, n_faces_of_nodes = calculate_face_properties.ordenate_faces_of_nodes_xy_plane(faces_centroids, faces_adj_by_nodes, nodes_centroids)
        
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
        self.data['n_faces_of_nodes'] = n_faces_of_nodes
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
          
    def create_3d_mesh_data(self) -> MeshProperty:
        self.init_mesh()
        self.init_3d_mesh_entities()
        self.create_3d_initial_array_properties()
        mesh_property = MeshProperty()
        mesh_property.insert_mesh_name(self.data['mesh_name'][0])
        mesh_property.insert_data(copy.deepcopy(self.data))
        mesh_property.export_data()
        return mesh_property
        
    def create_2d_mesh_data(self) -> MeshProperty:
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