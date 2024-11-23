# from ..running.run_simulation import RunSimulation
from packs.utils.info_manager import InfoManager
import os

def rodar():
    os.system('python3 testting1.py')

data_loaded = InfoManager('input_cards/inputs0_2.yml', 'input_cards/inputs0.yml')

data_loaded['gravity'] = False
data_loaded['load_mesh'] = False
data_loaded['convert_english_to_SI'] = False
data_loaded['n_test'] = 0
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = False
data_loaded['load_mesh'] = True
data_loaded['convert_english_to_SI'] = False
data_loaded['n_test'] = 1
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = False
data_loaded['load_mesh'] = False
data_loaded['convert_english_to_SI'] = True
data_loaded['n_test'] = 2
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = False
data_loaded['load_mesh'] = True
data_loaded['convert_english_to_SI'] = True
data_loaded['n_test'] = 3
data_loaded.save_obj()
rodar()

###################
data_loaded['gravity'] = True
data_loaded['load_mesh'] = False
data_loaded['convert_english_to_SI'] = False
data_loaded['n_test'] = 4
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = True
data_loaded['load_mesh'] = True
data_loaded['convert_english_to_SI'] = False
data_loaded['n_test'] = 5
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = True
data_loaded['load_mesh'] = False
data_loaded['convert_english_to_SI'] = True
data_loaded['n_test'] = 6
data_loaded.save_obj()
rodar()

data_loaded['gravity'] = True
data_loaded['load_mesh'] = True
data_loaded['convert_english_to_SI'] = True
data_loaded['n_test'] = 7
data_loaded.save_obj()
rodar()
