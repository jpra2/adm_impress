import numpy as np
from .impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh

class VugGenerator(object):
    def __init__(self, mesh_file, ellipsis_params_range, num_ellipsoids=10):
        self.mesh = FineScaleMesh(mesh_file)
        self.ellipsis_params_range = ellipsis_params_range
        self.num_ellipsoids = num_ellipsoids

    def run(self):
        # TODO: Add random rotation matrix.
        centroids = self.mesh.volumes.center[:]
        xs, ys, zs = centroids[:, 0], centroids[:, 1], centroids[:, 2]
        x_range = xs.min(), xs.max()
        y_range = ys.min(), ys.max()
        z_range = zs.min(), zs.max()
        centers, params = self.get_random_ellipsoids(x_range, y_range, z_range)

        # Compute vugs.
        # REVIEW: Maybe do this as vector/matrix operations.
        for center, param in zip(centers, params):
            vols_in_vug = ((xs - center[0]) / param[0])**2 \
                        + ((ys - center[1]) / param[1])**2 \
                        + ((zs - center[2]) / param[2])**2
            self.mesh.vug[vols_in_vug < 1] = 1

    def write_file(self, path="results/vugs.vtk"):
        vugs_meshset = self.mesh.core.mb.create_meshset()
        self.mesh.core.mb.add_entities(vugs_meshset, self.mesh.core.all_volumes)
        self.mesh.core.mb.write_file(path, [vugs_meshset])
    
    def get_random_ellipsoids(self, x_range, y_range, z_range):
        rng = np.random.default_rng()
        random_centers = np.zeros((self.num_ellipsoids, 3))
        random_params = np.zeros((self.num_ellipsoids, 3))

        random_centers[:, 0] = rng.uniform(low=x_range[0], high=x_range[1], size=self.num_ellipsoids)
        random_centers[:, 1] = rng.uniform(low=y_range[0], high=y_range[1], size=self.num_ellipsoids)
        random_centers[:, 2] = rng.uniform(low=z_range[0], high=z_range[1], size=self.num_ellipsoids)
        random_params[:] = rng.uniform(low=self.ellipsis_params_range[0], \
                                        high=self.ellipsis_params_range[1], \
                                        size=random_params.shape)
        
        return random_centers, random_params
