## criado por joao paulo
'''
classe Data criada para armazenar um dict linkando o nome da variavel
a variavel do impress e tambem as informacoes dos dados.
Essa classe será utilizada para automatizar alguns processos.
'''
import numpy as np
from . import directories_impress as direc
import pickle
import pdb

class Data:
    '''
    armazena as variaveis do impress e dicts que linkam
    as entidades do moab aos seus respectivos ids globais
    '''
    # valores_para_converter = ['hs', 'permeability', 'dist_cent']

    def __init__(self, n_nodes, n_edges, n_faces, n_volumes, fine_scale_mesh_obj):
        '''
        n_nodes: numero de nos
        n_edges: numero de arestas
        n_faces: numero de faces
        n_volumes: numero de volumes
        fine_scale_mesh_obj: objeto multiscalemeshMS
        '''

        self.variables = dict()
        self.info_data = dict()
        self.len_entities = dict()
        self.variables_impress = dict()

        self.mesh = fine_scale_mesh_obj
        self.name_variables = direc.path_local_variables
        self.name_info_data = direc.path_local_info_data
        fine_scale_mesh_obj.data = self

    def get_info_data(self, name, data_size, data_format, entity, level=0):
        '''
        name: nome da variavel
        impress_variable: variavel do impress
        data_size: tamanho da variavel
        data_format: formato da variavel
        entity: entidade onde foi setada
        level: nivel

        self.info_data: dict que armazena as informacoes das variaveis do impress
        '''
        # self.variables[name] = impress_variable
        info = dict()

        info[direc.names_datas[0]] = int(data_size)
        info[direc.names_datas[1]] = data_format
        info[direc.names_datas[2]] = entity
        info[direc.names_datas[3]] = level

        self.info_data[name] = info

    def init_datas(self, names=None):
        '''
        zera todas as variaveis do impress e armazena em self.variables
        names deve ser obrigatoriamente uma lista de strings
        '''

        variables = dict()
        variables_impress = dict()

        for name, infos in self.info_data.items():
            n = infos[direc.names_datas[0]]
            format = infos[direc.names_datas[1]]
            entity = infos[direc.names_datas[2]]
            n_entity = self.len_entities[entity]

            if format == direc.data_formats[0]:
                data = np.zeros(n_entity)
            elif format == direc.data_formats[1]:
                data = np.zeros(n_entity, dtype=np.int32)
            if n > 1:
                data = np.repeat(data, n).reshape([n_entity, n])

            variables[name] = data
            variables_impress[name] = name

        self.variables = variables
        self.variables_impress = variables_impress

    def init_dicts(self):
        '''
        dict_elements = dict com os dicts que linkam as entidades do pymoab aos
        seus ids globais
        mesma coisa para os centroides
        '''


        self.dict_elements = dict()

        try:
            dict_elem = dict(zip(self.mesh.core.all_volumes, self.mesh.volumes.all))
            self.dict_elements[direc.entities_lv0[3]] = dict_elem
        except:
            pass

        dict_elem = dict(zip(self.mesh.core.all_faces, self.mesh.faces.all))
        self.dict_elements[direc.entities_lv0[2]] = dict_elem
        dict_elem = dict(zip(self.mesh.core.all_edges, self.mesh.edges.all))
        self.dict_elements[direc.entities_lv0[1]] = dict_elem
        dict_elem = dict(zip(self.mesh.core.all_nodes, self.mesh.nodes.all))
        self.dict_elements[direc.entities_lv0[0]] = dict_elem

        self.elements_lv0 = dict()

        internal_faces = self.mesh.faces.internal
        self.elements_lv0[direc.entities_lv0_0[0]] = internal_faces
        vols_viz_faces = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 2, 3)
        self.elements_lv0[direc.entities_lv0_0[1]] = vols_viz_faces

        vols_viz_internal_faces = self.mesh.faces.bridge_adjacencies(internal_faces, 2, 3)
        self.elements_lv0[direc.entities_lv0_0[2]] = vols_viz_internal_faces

        boundary_faces = np.setdiff1d(self.mesh.faces.all, internal_faces)
        self.elements_lv0[direc.entities_lv0_0[4]] = boundary_faces

        vols_viz_boundary_faces = self.mesh.faces.bridge_adjacencies(boundary_faces, 2, 3).flatten()

        self.elements_lv0[direc.entities_lv0_0[3]] = vols_viz_boundary_faces

        self.elements_lv0[direc.entities_lv0[0]] = self.mesh.nodes.all
        self.elements_lv0[direc.entities_lv0[1]] = self.mesh.edges.all
        self.elements_lv0[direc.entities_lv0[2]] = self.mesh.faces.all
        self.elements_lv0[direc.entities_lv0[3]] = self.mesh.volumes.all

        # self.centroids = dict()
        # '''
        # self.centroids = dicts para linkar o nome das entidades aos seus centroides
        # '''
        #
        # try:
        #     self.centroids[direc.entities_lv0[3]] = self.mesh.volumes.center(self.mesh.volumes.all)
        # except:
        #     pass
        #
        # self.centroids[direc.entities_lv0[2]] = self.mesh.faces.center(self.mesh.faces.all)
        # self.centroids[direc.entities_lv0[1]] = self.mesh.edges.center(self.mesh.edges.all)
        # self.centroids[direc.entities_lv0[0]] = self.mesh.nodes.center(self.mesh.nodes.all)
        # self.variables[self.variables_impress['u_normal']] = np.absolute(self.mesh.faces.normal[:])
        # self.variables[self.variables_impress['NODES']] = self.centroids[direc.entities_lv0[0]].copy()
        # self.variables['u_normal'] = np.absolute(self.mesh.faces.normal[:])
        # self.variables['NODES'] = self.centroids[direc.entities_lv0[0]].copy()

    def update_variables_to_mesh(self, names=None):



        if names:
            for name in names:
                command = 'self.mesh.' + name + '[:] = ' + 'self.variables["' + name + '"]'
                exec(command)

        else:
            for name in self.variables.keys():
                command = 'self.mesh.' + name + '[:] = ' + 'self.variables["' + name + '"]'
                exec(command)

    def load_variables_from_mesh(self, names=None):

        if names:
            for name in names:
                command = 'self.variables["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

        else:
            for name in self.variables.keys():
                command = 'self.variables["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

    def export_variables_to_npz(self, file_name=None):

        name_variables = self.name_variables
        name_info_data = self.name_info_data

        if file_name:
            np.savez(file_name, **self.variables)
        else:
            np.savez(name_variables, **self.variables)

        # with open(self.name_info_data, 'rb') as f:
        #     pickle.dump

    def load_variables_from_npz(self, file_name=None):
        name_variables = self.name_variables
        name_info_data = self.name_info_data

        from .. import directories as direc

        if file_name:
            arq = np.load(file_name)
        else:
            arq = np.load(name_variables)

        for name, variable in arq.items():
            self.variables[name] = variable
            self.variables_impress[name] = name

        for key, value in direc.variables_impress.items():
            self.variables_impress[key] = value

    def save_info_data(self):

        name_info_data = self.name_info_data

        with open(name_info_data, 'wb') as f:
            pickle.dump(self.info_data, f)

    def load_info_data(self):

        name_info_data = self.name_info_data

        with open(name_info_data, 'rb') as f:
            self.info_data = pickle.loads(f.read())

    def __setitem__(self, key, value):
        self.variables[self.variables_impress[key]] = value

    def __getitem__(self, key):
        return self.variables[self.variables_impress[key]]

    def __hash__(self, key):
        return hash(self.variables[self.variables_impress[key]])
