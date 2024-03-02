from packs.manager.arraydatamanager import SuperArrayManager
import meshio
import numpy as np

class MeshioWrapper:

    def __init__(self, mesh_path):
        self.msh: meshio._mesh.Mesh = meshio.read(mesh_path)

    @property
    def points_centroids(self):
        return self.msh.points
    
    @property
    def points(self):
        n_points = self.points_centroids.shape[0]
        return np.arange(n_points)
    
    @property
    def lines_points(self):
        data = self.msh.cells_dict
        return data['line']
    
    @property
    def lines(self):
        n_lines = self.lines_points.shape[0]
        return np.arange(n_lines)
    
    @property
    def triangles_points(self):
        data = self.msh.cells_dict
        return data['triangle']
    
    @property
    def triangles(self):
        n_tris = self.triangles_points.shape[0]
        return np.arange(n_tris)
    
    @property
    def quads_points(self):
        data = self.msh.cells_dict
        return data['quad']
    
    @property
    def quads(self):
        n_quads = self.quads_points.shape[0]
        return np.arange(n_quads)

        
