import yaml
import os
import numpy as np
from pymoab import core, types, rng, topo_util

class createMesh:

    def __init__(self, dim=3):
        input_file = os.path.join('input_cards', 'generate_structured_mesh.yml')

        with open(input_file, 'r') as f:
            self.mesh_data = yaml.safe_load(f)

        self.mesh_data['block_size'] = np.array(self.mesh_data['block_size'])
        self.mesh_data['mesh_size'] = np.array(self.mesh_data['block_size'])*np.array(self.mesh_data['block_number'])
        self.mesh_data['starting_point'] = np.array(self.mesh_data['starting_point'])
        self.dim = int(self.mesh_data['dimension'])
        assert self.dim in [2, 3]

        self.mb = core.Core()
        self.root_set = self.mb.get_root_set()
        self.mtu = topo_util.MeshTopoUtil(self.mb)

    def init_params(self):
        
        block_size = self.mesh_data['block_size']
        mesh_size = self.mesh_data['mesh_size']

        self.params = dict()
        # self.params['nblocks'] = np.floor(mesh_size/block_size).astype(int)
        self.params['nblocks'] = np.array(self.mesh_data['block_number'])
        nblocks = self.params['nblocks']

        for i in range(3):
            test = nblocks[i]*block_size[i]
            if test != mesh_size[i]:
                raise ValueError('Dados de malha nao servem para malha estruturada')

    def create_fine_vertices(self):
        if self.dim == 3:
            self.create_fine_vertices_3D()
        elif self.dim == 2:
            self.create_fine_vertices_2D()


    def create_fine_vertices_3D(self):

        coords = np.array([(i, j, k)
                           for k in (
                               np.arange(
                                   self.params['nblocks'][2]+1, dtype='float64') *self.mesh_data['block_size'][2])
                           for j in (
                               np.arange(
                                   self.params['nblocks'][1]+1, dtype='float64') *self.mesh_data['block_size'][1])
                           for i in (
                               np.arange(
                                   self.params['nblocks'][0]+1, dtype='float64') *self.mesh_data['block_size'][0])
                           ], dtype='float64')
        coords+=self.mesh_data['starting_point']
        self.verts = self.mb.create_vertices(coords.flatten())

    def _create_hexa(self, i, j, k):
    
        hexa = [self.verts[int((i)+(j*(self.params['nblocks'][0]+1))+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i, j, k)
                self.verts[int((i+1)+(j*(self.params['nblocks'][0]+1))+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i+1, j, k)
                self.verts[int((i+1)+(j+1)*(self.params['nblocks'][0])+(j+1)+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i+1, j+1, k)
                self.verts[int((i)+(j+1)*(self.params['nblocks'][0])+(j+1)+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i, j+1, k)

                self.verts[int((i)+(j*(self.params['nblocks'][0]+1))+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i, j, k+1)
                self.verts[int((i+1)+(j*(self.params['nblocks'][0]+1))+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i+1, j, k+1)
                self.verts[int((i+1)+(j+1)*(self.params['nblocks'][0])+(j+1)+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))],  # (i+1, j+1, k+1)
                self.verts[int((i)+(j+1)*(self.params['nblocks'][0])+(j+1)+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1))))]]  # (i, j+1, k+1)
        
        return hexa

    def _create_quads(self, i, j):
        quad = [self.verts[int((i)+(j*(self.params['nblocks'][0]+1)))],  # (i, j)
                self.verts[int((i+1)+(j*(self.params['nblocks'][0]+1)))],  # (i+1, j)
                self.verts[int((i+1)+(j+1)*(self.params['nblocks'][0])+(j+1))],  # (i+1, j+1)
                self.verts[int((i)+(j+1)*(self.params['nblocks'][0])+(j+1))] # (i, j+1)
          ]

        return quad  

    def create_elements(self):
        if self.dim == 3:
            self.create_elements_3D()
        elif self.dim == 2:
            self.create_elements_2D()

    def create_elements_3D(self):
        nbs = self.params['nblocks']
        hexas = [self._create_hexa(i, j, k) for i in range(nbs[0]) for j in range(nbs[1]) for k in range(nbs[2])]
        self.mb.create_elements(types.MBHEX, hexas)

    def export_mesh(self):
        mesh_name = os.path.join('mesh', self.mesh_data['mesh_name'])
        self.mb.write_file(mesh_name)
        self.mb.write_file(os.path.join('mesh', 'test_mesh.vtk'))

    def create_fine_vertices_2D(self):
        coords = np.array([(i, j, 0)
                           for j in (
                               np.arange(
                                   self.params['nblocks'][1]+1, dtype='float64') *self.mesh_data['block_size'][1])
                           for i in (
                               np.arange(
                                   self.params['nblocks'][0]+1, dtype='float64') *self.mesh_data['block_size'][0])
                           ], dtype='float64')
        coords+=self.mesh_data['starting_point']
        self.verts = self.mb.create_vertices(coords.flatten())
    
    def create_elements_2D(self):
        nbs = self.params['nblocks'][0:2]
        quads = [self._create_quads(i, j) for i in range(nbs[0]) for j in range(nbs[1])]
        self.mb.create_elements(types.MBQUAD, quads)

    def generate_non_uniform_mesh(self):
        if self.dim == 2:
            self.generate_non_uniform_mesh_2D()
        else:
            raise ValueError
    

    def generate_non_uniform_mesh_2D(self):
        # x_points = np.array(self.mesh_data['nonuniform_mesh']['x_points']).astype(np.float64)
        # y_points = np.array(self.mesh_data['nonuniform_mesh']['y_points'])

        x_data = self.mesh_data['nonuniform_mesh']['x_points']
        y_data = self.mesh_data['nonuniform_mesh']['y_points']

        x_points = np.array([eval(i) for i in x_data])
        y_points = np.array([eval(i) for i in y_data])

        self.params = dict()
        nbs = np.array([x_points.shape[0] - 1, y_points.shape[0] - 1, 0])
        self.params['nblocks'] = nbs

        coords = np.array([(i, j, 0)
                           for j in y_points
                           for i in x_points
                           ], dtype=np.float64)
        
        self.verts = self.mb.create_vertices(coords.flatten())
        quads = [self._create_quads(i, j) for i in range(nbs[0]) for j in range(nbs[1])]
        self.mb.create_elements(types.MBQUAD, quads)