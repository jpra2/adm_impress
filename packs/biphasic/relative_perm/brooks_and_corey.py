import numpy as np
from typing import Tuple

class BrooksAndCorey:

    def __init__(self, Sor: float= 0.2, Swc: float=0.2, nw: int=2, no: int=2, krw0: float=1, kro0: float=1):
        """Biphasic Corey and Brooks relative permeability Init parameters

        Args:
            Sor (float): residual oil saturation
            Swc (float): critical water saturation
            nw (int): water exponent
            no (int): oil exponent
            krw0 (float): initial water relative permeability
            kro0 (float): initial oil relative permeability
        """

        self.Sor = Sor
        self.Swc = Swc
        self.n_w = nw
        self.n_o = no
        self.krw0 = krw0
        self.kro0 = kro0

    def _stemp(self, S:np.ndarray) -> np.ndarray:
        """Update self.stemp

        Args:
            S (np.ndarray): water saturation

        Returns:
            (np.ndarray): stemp
        """
        S1 = S.copy()

        S1[S>1 - self.Sor] = 1 - self.Sor
        S1[S<self.Swc] = self.Swc
        return (S1 - self.Swc) / (1 - self.Swc - self.Sor)

    def _krw(self, S_temp: np.ndarray) -> np.ndarray:
        return self.krw0*(np.power(S_temp, self.n_w))

    def _kro(self, S_temp: np.ndarray) -> np.ndarray:
        return self.kro0*(np.power(1 - S_temp, self.n_o))

    def calculate(self, saturations: np.ndarray, test_saturation=True) -> Tuple[np.ndarray, np.ndarray]:
        
        if test_saturation is True:
            self._test_saturations(saturations)
        stemp = self._stemp(saturations)
        krw = self._krw(stemp)
        kro = self._kro(stemp)

        return krw, kro

    def _test_saturations(self, saturations: np.ndarray) -> None:
        """Test if the saturations is between max and min values 

        Args:
            saturations (np.ndarray): water saturations
        """
        # test1 = saturations < self.Swc
        # test2 = saturations > 1 - self.Sor
        test1 = saturations < 0
        test2 = saturations > 1
        n = test1.sum() + test2.sum()
        assert n == 0

    def __call__(self, saturations:np.ndarray) -> np.ndarray:
        return self.calculate(saturations)
