import numpy as np
from typing import Tuple
from packs.manager import BoundaryConditions

class BiphasicMobility:

    def __init__(self, miw: float=1.0, mio: float=1.2):
        """Initialize parameters

        Args:
            miw (float): water dynamic viscosity
            mio (float): oil dynamic viscosity
        """
        self.miw = miw
        self.mio = mio

    def mobw(self, krw):
        return krw/self.miw

    def mobo(self, kro):
        return kro/self.mio        
    
    def calculate(self, krw: np.ndarray, kro: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """_summary_

        Args:
            krw (np.ndarray): water relative permeability
            kro (np.ndarray): oil relative permeability

        Returns:
            tuple: water and oil mobility
        """

        # lambda_w = krw/self.miw
        # lambda_o = kro/self.mio
        return self.mobw(krw), self.mobo(kro)
    
    def get_total_mobility(self, water_mobility: np.ndarray, oil_mobility: np.ndarray) -> np.ndarray:
        return water_mobility + oil_mobility

    def update_edges_saturation_foum(self, faces_saturation:np.ndarray, edges_flux: np.ndarray, adjacencies: np.ndarray, bc: BoundaryConditions, edges_saturation: np.ndarray, bool_boundary_edges: np.ndarray) -> np.ndarray:
        
        test = edges_flux >= 0
        edges_saturation[test] = faces_saturation[adjacencies[test, 0]]
        test[:] = ~test
        edges_saturation[test] = faces_saturation[adjacencies[test, 1]]
        edges_saturation[bool_boundary_edges] = faces_saturation[adjacencies[bool_boundary_edges, 0]]

        edges_sat_presc = bc['water_saturation_edges']['id']
        if len(edges_sat_presc) > 0:
            sat_presc_value = bc['water_saturation_edges']['value']
            edges_saturation[edges_sat_presc] = sat_presc_value
        
        return edges_saturation

    def get_fw(self, water_mobility, oil_mobility):
        total_mobility = self.get_total_mobility(water_mobility, oil_mobility)
        return water_mobility/total_mobility




        


