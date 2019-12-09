# from ..running.run_simulation import RunSimulation
from packs.utils.info_manager import InfoManager
import os

data_loaded = InfoManager('input_cards/inputs0.yml')

data_loaded['gravity'] = False
data_loaded['load_mesh'] = False
data_loaded['convert_english_to_SI'] = False
data_loaded['n_test'] = 0
data_loaded.save_obj()

os.system('python3 ')
