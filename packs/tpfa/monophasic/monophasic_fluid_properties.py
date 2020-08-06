from packs.data_class import DataManager
import yaml

class MonophasicFluidProperties:
    __slots__ = ['_density', '_mi']

    def __init__(self, input_file=''):

        if input_file == '':
            input_file = 'input_cards/monophasic/monophasic_fluid_properties.yml'

        with open(input_file, 'r') as f:
            data_loaded = yaml.safe_load(f)

        self._density = float(data_loaded['monophasic_data']['density'])
        self._mi = float(data_loaded['monophasic_data']['mi'])

    @property
    def density(self):
        return self._density

    @property
    def mi(self):
        return self._mi
