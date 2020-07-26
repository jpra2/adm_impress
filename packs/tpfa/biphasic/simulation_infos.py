import os
import yaml
from copy import deepcopy
import numpy as np

class SimulationInfos:

    def __init__(self):
        path = os.path.join('input_cards', 'biphasic')
        path = os.path.join(path, 'simulation_info.yml')

        with open(path, 'r') as f:
            data_loaded = yaml.safe_load(f)

        for name in data_loaded.keys():
            try:
                data_loaded[name] = float(data_loaded[name])
            except:
                pass
        data_loaded['show_vtk_for_vpis'] = np.array(data_loaded['show_vtk_for_vpis'], dtype=float)

        self.data = deepcopy(data_loaded)
