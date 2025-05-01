import meshio
import numpy as np
from packs.utils.test_functions import test_mesh_path
from copy import deepcopy

class MeshioWrapper:
    tags_to_remove = ['gmsh:bounding_entities']

    def __init__(self, mesh_path):   
        full_path = test_mesh_path(mesh_path)
        self.msh: meshio._mesh.Mesh = meshio.read(full_path)

    @property
    def points_centroids(self) -> np.ndarray:
        return self.msh.points
    
    @property
    def points(self) -> np.ndarray:
        n_points = self.points_centroids.shape[0]
        return np.arange(n_points)
    
    @property
    def lines_points(self) -> np.ndarray:
        data = self.msh.cells_dict
        return data['line']
    
    @property
    def lines_centroids(self) -> np.ndarray:
        centroids_lines_points = self.points_centroids[self.lines_points]
        l_centroids = np.mean(centroids_lines_points, axis=1)
        return l_centroids
    
    @property
    def lines(self) -> np.ndarray:
        n_lines = self.lines_points.shape[0]
        return np.arange(n_lines)
    
    @property
    def triangles_points(self) -> np.ndarray:
        data = self.msh.cells_dict
        return data['triangle']
    
    @property
    def triangles(self) -> np.ndarray:
        n_tris = self.triangles_points.shape[0]
        return np.arange(n_tris)
    
    @property
    def quads_points(self) -> np.ndarray:
        data = self.msh.cells_dict
        return data['quad']
    
    @property
    def quads(self) -> np.ndarray:
        n_quads = self.quads_points.shape[0]
        return np.arange(n_quads)

    @property
    def physical_tags(self) -> list:
        tags = set(list(self.msh.cell_sets_dict.keys())) - set(self.tags_to_remove)
        return list(tags)

    @property
    def physical_int_tags(self):
        try:
            return deepcopy(self.msh.cell_data_dict["gmsh:physical"])
        except KeyError:
            return dict()

    def get_elements_by_physical_tag(self, tag: str) -> dict:
        return deepcopy(self.msh.cell_sets_dict[tag])
    
    def get_elements_by_physical_int_tag(self, key: str, tag: int):

        data = self.physical_int_tags
        test = data[key] == tag
        
        if key == 'vertex':
            vertices = self.msh.cells_dict[key].flatten()
            return vertices[test]
        elif key == 'line':
            lines = self.lines
            return lines[test]
        elif key == 'quad':
            return self.quads[test]
        elif key == 'triangle':
            return self.triangles[test]
        else:
            raise ValueError
        
