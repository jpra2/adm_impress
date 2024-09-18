import numpy as np
from typing import Tuple

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

    def get_water_fractionary_flux(self, water_mobility: np.ndarray, oil_mobility: np.ndarray) -> np.ndarray:
        return water_mobility/self.get_total_mobility(water_mobility, oil_mobility)