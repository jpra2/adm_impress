import yaml
import os
import numpy as np
from pymoab import core, types, rng, topo_util


class createNonUniformMesh:

    def __init__(self):
        self.mb = core.Core()
        self.root_set = self.mb.get_root_set()
        self.mtu = topo_util.MeshTopoUtil(self.mb)

    def create_fine_vertices(self, x_points, y_points, z_points):
        npoints = np.array([len(x_points), len(y_points), len(z_points)], dtype=int)

        coords = np.array([(i, j, k)
                           for k in (z_points)
                           for j in (y_points)
                           for i in (x_points)
                           ], dtype='float64')

        verts = self.mb.create_vertices(coords.flatten())
        return verts, npoints

    def _create_hexa(self, i, j, k, verts, nblocks):
        hexa = [verts[(i)+(j*(nblocks[0]+1))+(k*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i, j, k)
                verts[(i+1)+(j*(nblocks[0]+1))+(k*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i+1, j, k)
                verts[(i+1)+(j+1)*(nblocks[0])+(j+1)+(k*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i+1, j+1, k)
                verts[(i)+(j+1)*(nblocks[0])+(j+1)+(k*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i, j+1, k)

                verts[(i)+(j*(nblocks[0]+1))+((k+1)*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i, j, k+1)
                verts[(i+1)+(j*(nblocks[0]+1))+((k+1)*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i+1, j, k+1)
                verts[(i+1)+(j+1)*(nblocks[0])+(j+1)+((k+1)*((nblocks[0]+1)*(nblocks[1]+1)))],  # (i+1, j+1, k+1)
                verts[(i)+(j+1)*(nblocks[0])+(j+1)+((k+1)*((nblocks[0]+1)*(nblocks[1]+1)))]]  # (i, j+1, k+1)

        return hexa

    def create_elements(self, verts, npoints):
        nblocks = npoints-1
        hexas = [self._create_hexa(i, j, k, verts, nblocks) for i in range(nblocks[0]) for j in range(nblocks[1]) for k in range(nblocks[2])]
        self.mb.create_elements(types.MBHEX, hexas)

    def export_mesh(self, mesh_name):
        self.mb.write_file(mesh_name)

    def create_mesh(self, x_points, y_points, z_points, mesh_name):
        verts, npoints = self.create_fine_vertices(x_points, y_points, z_points)
        self.create_elements(verts, npoints)
        self.export_mesh(mesh_name)

        print(f'\n{mesh_name} created\n')
