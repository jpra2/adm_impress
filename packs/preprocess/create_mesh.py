import yaml
import os
import numpy as np
from pymoab import core, types, rng, topo_util

class createMesh:

    def __init__(self):
        input_file = os.path.join('input_cards', 'generate_structured_mesh.yml')

        with open(input_file, 'r') as f:
            self.mesh_data = yaml.safe_load(f)

        self.mesh_data['block_size'] = np.array(self.mesh_data['block_size'])
        self.mesh_data['mesh_size'] = np.array(self.mesh_data['block_size'])*np.array(self.mesh_data['block_number'])
        self.mesh_data['starting_point'] = np.array(self.mesh_data['starting_point'])
        self.init_params()

        self.mb = core.Core()
        self.root_set = self.mb.get_root_set()
        self.mtu = topo_util.MeshTopoUtil(self.mb)

    def init_params(self):
        starting_point = self.mesh_data['starting_point']
        block_size = self.mesh_data['block_size']
        mesh_size = self.mesh_data['mesh_size']


        self.params = dict()
        self.params['nblocks'] = np.floor(mesh_size/block_size).astype(int)
        nblocks = self.params['nblocks']

        for i in range(3):
            test = nblocks[i]*block_size[i]
            if test != mesh_size[i]:
                raise ValueError('Dados de malha nao servem para malha estruturada')

    def create_fine_vertices(self):

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
        hexa = [self.verts[(i)+(j*(self.params['nblocks'][0]+1))+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i, j, k)
                self.verts[(i+1)+(j*(self.params['nblocks'][0]+1))+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i+1, j, k)
                self.verts[(i+1)+(j+1)*(self.params['nblocks'][0])+(j+1)+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i+1, j+1, k)
                self.verts[(i)+(j+1)*(self.params['nblocks'][0])+(j+1)+(k*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i, j+1, k)

                self.verts[(i)+(j*(self.params['nblocks'][0]+1))+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i, j, k+1)
                self.verts[(i+1)+(j*(self.params['nblocks'][0]+1))+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i+1, j, k+1)
                self.verts[(i+1)+(j+1)*(self.params['nblocks'][0])+(j+1)+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))],  # (i+1, j+1, k+1)
                self.verts[(i)+(j+1)*(self.params['nblocks'][0])+(j+1)+((k+1)*((self.params['nblocks'][0]+1)*(self.params['nblocks'][1]+1)))]]  # (i, j+1, k+1)

        return hexa

    def create_elements(self):
        nbs = self.params['nblocks']
        hexas = [self._create_hexa(i, j, k) for i in range(nbs[0]) for j in range(nbs[1]) for k in range(nbs[2])]
        self.mb.create_elements(types.MBHEX, hexas)

    def export_mesh(self):
        mesh_name = os.path.join('mesh', self.mesh_data['mesh_name'])
        self.mb.write_file(mesh_name)
