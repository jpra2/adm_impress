from packs.tpfa.biphasic import relative_permeability
import yaml
import pdb

def get_relative_permeability(input_file=''):

    if input_file == '':
        input_file = 'input_cards/biphasic/relative_permeability_data.yml'

    with open(input_file, 'r') as f:
        data_loaded = yaml.safe_load(f)

    name_data = data_loaded['relative_permeability']
    obj = getattr(relative_permeability, name_data)
    obj = obj(data_loaded)
    return obj
