from packs.utils.utils_old import get_box
from packs.data_class.data_manager import DataManager
from packs.directories import only_mesh_name
from packs.properties import Data
import numpy as np

class RockProperties:

    _type_set_permeability = ['all']
    _type_permeability = ['symetric']
    _n_volumes_kwarg = ['n_volumes']
    _dimensions = [1, 2, 3]
    _properties = ['porosity', 'permeability']
    _perm_loc = ['kxx', 'kxy', 'kxz', 'kyx', 'kyy', 'kyz', 'kzx', 'kzy', 'kzz']
    _kwargs = _type_set_permeability + _type_permeability + _n_volumes_kwarg + _dimensions + _properties + _perm_loc

    campos_set_permeability = ['type_set_permeability', 'indices', 'value']

    '''
        The RockProperties class
    '''

    def __init__(self, load=False, dimension=3, n_volumes=0):
        data_name = type(self).__name__
        data_name = data_name + '_' + only_mesh_name + '.npz'
        self.dimension = dimension
        if load == True:
            self.n_volumes = len(self['permeability'])
        else:
            self.n_volumes = n_volumes
            self.init_permeability_matrix()

    def set_permeability(self, *args, **kwargs):
        campos1 = []
        for i, name in enumerate(RockProperties.campos_set_permeability):
            campos1.append(kwargs.get(name))
        campos2 = []
        for i, name in enumerate(RockProperties._perm_loc):
            campos2.append(kwargs.get(name))

        self._set_permeability(campos1, campos2)

    def _set_permeability(self, campos1, campos2):
        campos2_ = [i if i != None else 0.0 for i in campos2]

        if campos1[0] == 'all':
            permeability = np.array(campos2_).reshape((self.dimension, self.dimension))
            self['permeability'] = np.repeat(permeability, self.n_volumes)
        else:
            indices = campos1[1]
            values = campos1[2]
            self['permeability'][indices] = values

    def set_porosity(self, porosity):
        self['porosity'] = porosity

    def init_permeability_matrix(self):
        self.permeability = Data('permeability')
        self.permeability = np.zeros((self.n_volumes, self.dimension, self.dimension))

    @property
    def permeability(self):
        try:
            return self['permeability']
        except AttributeError:
            print('\nPermeability not found\n')
