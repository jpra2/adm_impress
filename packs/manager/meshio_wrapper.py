import meshio
import numpy as np

class MeshioWrapper:

    def __init__(self, mesh_path):
        self.msh: meshio._mesh.Mesh = meshio.read(mesh_path)

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
        tags = list(self.msh.cell_sets_dict.keys())
        return tags

    def get_elements_by_physical_tag(self, tag: str) -> dict:
        return self.msh.cell_sets_dict[tag]



        
