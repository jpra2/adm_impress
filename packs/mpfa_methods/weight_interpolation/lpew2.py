import numpy as np
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation

class Lpew2Weight:

    @staticmethod
    def preprocess(
        nodes_centroids,
        unitary_normal_edges,
        nodes_of_edges,
        edges, 
        **kwargs
    ):

        return LsdsFluxCalculation.define_A_B_points_of_edges(
            nodes_centroids,
            unitary_normal_edges,
            nodes_of_edges,
            edges
        )
    
    

