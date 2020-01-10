## criado por joao paulo
'''
classe Data criada para armazenar um dict linkando o nome da variavel
a variavel do impress e tambem as informacoes dos dados.
Essa classe será utilizada para automatizar alguns processos.
'''
import numpy as np
from .. import directories as direc
import pickle
import pdb
from .data_manager import DataManager
from ..convert_unit import constants

class Data(DataManager):
    '''
    armazena as variaveis do impress e dicts que linkam
    as entidades do moab aos seus respectivos ids globais
    '''
    # valores_para_converter = ['hs', 'permeability', 'dist_cent']

    def __init__(self, fine_scale_mesh_obj, elementsLv0_obj, load: bool=False, data_name: str='data_impress.npz'):
        '''
        fine_scale_mesh_obj: objeto multiscalemeshMS
        '''
        super().__init__(data_name, load=load)

        self.info_data = dict()
        self.len_entities = dict()
        self.len_entities['faces'] = len(elementsLv0_obj['faces'])
        self.len_entities['nodes'] = len(elementsLv0_obj['nodes'])
        self.len_entities['volumes'] = len(elementsLv0_obj['volumes'])
        self.len_entities['edges'] = len(elementsLv0_obj['edges'])
        self.variables_impress = dict()

        self.mesh = fine_scale_mesh_obj
        fine_scale_mesh_obj.data = self
        if not load:
            self.run()

        if load:
            self.get_variables_impress()

        self._loaded = True

    def get_info_data(self):
        '''
        self.info_data: dict que armazena as informacoes das variaveis do impress
        '''
        # self.variables[name] = impress_variable
        for name, data in direc.variables_loaded.items():

            info = dict()
            info['data_size'] = data['data size']
            info['data_format'] = data['data type']
            info['entity'] = data['type']
            info['level'] = data['level']
            self.info_data[name] = info

    def init_datas(self):
        '''
        zera todas as variaveis do impress e armazena
        '''
        variables_impress = dict()

        for name, infos in self.info_data.items():
            n = infos['data_size']
            format = infos['data_format']
            entity = infos['entity']
            n_entity = self.len_entities[entity]

            if format == 'float':
                data = np.zeros(n_entity)
            elif format == 'int':
                data = np.zeros(n_entity, dtype=np.int32)
            if n > 1:
                data = np.repeat(data, n).reshape([n_entity, n])

            self[name] = data
            self._data[name] = data
            variables_impress[name] = name

        self.variables_impress = variables_impress

    def update_variables_to_mesh(self, names=None):

        if names:
            for name in names:
                command = 'self.mesh.' + name + '[:] = ' + 'self._data["' + name + '"]'
                exec(command)

        else:
            for name in self._data.keys():
                command = 'self.mesh.' + name + '[:] = ' + 'self._data["' + name + '"]'
                try:
                    exec(command)
                except:
                    print(command)
                    import pdb; pdb.set_trace()

    def load_variables_from_mesh(self, names=None):

        if names:
            for name in names:
                command = 'self._data["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

        else:
            for name in self._data.keys():
                command = 'self._data["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

    def get_variables_impress(self):
        self.variables_impress = dict()
        for name in self._data.keys():
            self.variables_impress[name] = name

    def run(self):
        self.get_info_data()
        self.init_datas()
