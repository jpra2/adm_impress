import numpy as np
from .. import directories as direc
import pdb
from ..directories import data_loaded
from .data_manager import DataManager
from ..convert_unit import constants

class FLuidProperties:
    def __init__(self, load: bool=False, data_name: str='fluid_data.npz'):
        compositional_info = data_loaded['compositional_data']
        self.z = compositional_info['z']
        self.Vc = compositional_info['Vc']
        self.Pc = compositional_info['Pc']
        self.Tc = compositional_info['Tc']
        self.Mw = compositional_info['Mw']
        self.P = compositional_info['P']
        self.T = compositional_info['T']
        self.w = compositional_info['w']
        self.Bin = compositional_info['Bin']
        self.Nc = len(Mw)
