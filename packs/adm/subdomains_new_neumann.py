import numpy as np

class SubdomainsNeumman:

    _attrs = ['volumes', 'pressure_volumes', 'pressure_values',
              'flux_volumes', 'flux_values', 'faces', 'internal', 'keq_internal_faces'
              'total_velocity_source_term_faces', 'total_flux_source_term_faces',
              'volumes_adj_internal_faces', 'volumes_adj_boundary_faces']

    def __init__(self):
        self.subdomains = []

    def reset(self):
        self.subdomains = []

    def insert(self, obj):
        assert list(obj.keys()) = SubdomainsNeumman._attrs
        self.subdomains.append(obj)

    def calculate_pressure(self):
