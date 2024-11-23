import yaml

class BiphasicFluidProperties:

    def __init__(self, input_file=''):

        if input_file == '':
            input_file = 'input_cards/biphasic/biphasic_fluid_properties.yml'

        with open(input_file, 'r') as f:
            data_loaded = yaml.safe_load(f)

        self._density_w = float(data_loaded['density_w'])
        self._density_o = float(data_loaded['density_o'])
        self._mi_w = float(data_loaded['mi_w'])
        self._mi_o = float(data_loaded['mi_o'])

    @property
    def rho_w(self):
        return self._density_w

    @property
    def rho_o(self):
        return self._density_o

    @property
    def mi_o(self):
        return self._mi_o

    @property
    def mi_w(self):
        return self._mi_w
