
class SimulationVariables:

    def gradient_pressure(self, pressure_volumes, centroid_volumes, volumes_adj_internal_faces):

        delta_p = self.delta_p(pressure_volumes, volumes_adj_internal_faces)
        direction = self.direction(centroid_volumes, volumes_adj_internal_faces)
        norm_direction = np.linalg.norm(direction, axis=1)
        direction = direction/norm_direction
        gradient_pressure = -(delta_p/norm_direction) * direction
        return gradient_pressure

    def delta_p(self, pressure_volumes, volumes_adj_internal_faces):
        delta_p = pressure_volumes[volumes_adj_internal_faces]
        return delta_p[:, 1] - delta_p[:, 0]

    def direction(self, centroid_volumes, volumes_adj_internal_faces):
        direction = centroid_volumes[volumes_adj_internal_faces]
        return direction[:, 1] - direction[:, 0]

    def total_velocity(self, pressure_volumes, centroid_volumes, volumes_adj_internal_faces, permeability, transmissibility_internal_faces):
        pass
