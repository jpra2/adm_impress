import numpy as np
from .impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh

class VugGenerator(object):
    def __init__(self, mesh_file, ellipsis_params_range, num_ellipsoids=10):
        """
        Constructor method.

        Inputs
        ------
        mesh_file (String): A string containing the path to the input mesh file.
        ellipsis_params_range (iterable of size 2): A iterable containing the maximum and the
        minimum value for the ellipsoids parameters.
        num_ellipsoids (int): the number of ellipsoids to create.

        """
        self.mesh = FineScaleMesh(mesh_file)
        self.ellipsis_params_range = ellipsis_params_range
        self.num_ellipsoids = num_ellipsoids

    def run(self):
        """
        Main method. Creates random sized and rotated ellipsoids and
        assigns the volumes inside the vugs.
        
        Inputs
        ------
        None

        Output
        ------
        None

        """
        centroids = self.mesh.volumes.center[:]
        xs, ys, zs = centroids[:, 0], centroids[:, 1], centroids[:, 2]
        x_range = xs.min(), xs.max()
        y_range = ys.min(), ys.max()
        z_range = zs.min(), zs.max()
        centers, params, angles = self.get_random_ellipsoids(x_range, y_range, z_range)

        # Compute vugs.
        for center, param, angle in zip(centers, params, angles):
            R = self.get_rotation_matrix(angle)
            X = (centroids - center).dot(R.T)
            vols_in_vug = (X / param)**2
            vols_in_vug = vols_in_vug.sum(axis=1)
            self.mesh.vug[vols_in_vug < 1] = 1

    def write_file(self, path="results/vugs.vtk"):
        """
        Writes the resultant mesh into a file. Default path is 'results/vugs.vtk'.

        Inputs
        ------
        path (String): A string containing the output file name.

        Outputs
        -------
        None
        
        """
        vugs_meshset = self.mesh.core.mb.create_meshset()
        self.mesh.core.mb.add_entities(vugs_meshset, self.mesh.core.all_volumes)
        self.mesh.core.mb.write_file(path, [vugs_meshset])
    
    def get_random_ellipsoids(self, x_range, y_range, z_range):
        """
        Generates random points as the ellipsoids centers as the axis sizes 
        and random rotation angles with respect to the cartesian coordinates (x,y,z).

        Inputs
        ------
        x_range (iterable of size 2): A iterable containing the maximum and minimum 
        values of the x coordinate.
        y_range (iterable of size 2): A iterable containing the maximum and minimum 
        values of the y coordinate.
        z_range (iterable of size 2): A iterable containing the maximum and minimum 
        values of the z coordinate.

        Outputs
        -------
        random_centers (num_ellipsoids x 3 numpy array): The generated center 
        points for the ellipsoids.
        random_params (num_ellipsoids x 3 numpy array): The parameters a.k.a 
        the size of the three axis of the ellipsoids.
        random_angles (num_ellipsoids x 3 numpy array): The rotation angles 
        for each ellipsoid.
        
        """
        rng = np.random.default_rng()
        random_centers = np.zeros((self.num_ellipsoids, 3))

        random_centers[:, 0] = rng.uniform(low=x_range[0], high=x_range[1], size=self.num_ellipsoids)
        random_centers[:, 1] = rng.uniform(low=y_range[0], high=y_range[1], size=self.num_ellipsoids)
        random_centers[:, 2] = rng.uniform(low=z_range[0], high=z_range[1], size=self.num_ellipsoids)
        random_params = rng.uniform(low=self.ellipsis_params_range[0], \
                                        high=self.ellipsis_params_range[1], \
                                        size=(self.num_ellipsoids, 3))
        random_angles = rng.uniform(low=0.0, high=2*np.pi, size=(self.num_ellipsoids, 3))
        
        return random_centers, random_params, random_angles
    
    def get_rotation_matrix(self, angle):
        """
        Calculates the rotation matrix for the given angle. Such matrix is the result of 
        three rotations: first with respect to x, then y, then z.

        Inputs
        ------
        angle (3 x 1 numpy array): The three rotation angles in radians.
        
        Outputs
        -------
        R (3 x 3 numpy array): The rotation matrix.

        """
        cos_ang = np.cos(angle)
        sin_ang = np.sin(angle)
        R = np.array([
            cos_ang[1]*cos_ang[2], -cos_ang[1]*sin_ang[2], sin_ang[1], \
            sin_ang[0]*sin_ang[1]*cos_ang[2] + cos_ang[0]*sin_ang[2], \
            cos_ang[0]*cos_ang[2] - sin_ang[0]*sin_ang[1]*sin_ang[2], \
            -sin_ang[0]*sin_ang[1], sin_ang[0]*sin_ang[2] - cos_ang[0]*sin_ang[1]*sin_ang[2], \
            cos_ang[0]*sin_ang[1]*sin_ang[2] + sin_ang[0]*cos_ang[2], cos_ang[0]*cos_ang[1]
        ]).reshape((3,3))

        return R
