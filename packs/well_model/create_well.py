import weakref
import os
from packs.directories import only_mesh_name
from packs.well_model import peaceman
import pandas as pd
import numpy as np


def test_id_text(id_text):
    if id_text in AllWells.wells.keys():
        raise NameError(f'id_text not must be in {AllWells.wells.keys()}')


class AllWells:
    wells = dict()
    database_name = os.path.join('flying', 'wells_database_' + only_mesh_name + '.csv')

    # database_name = os.path.join('flying', 'wells_database_' + only_mesh_name + '.h5')

    @classmethod
    def create_database(cls):
        list_ids = []
        list_volumes_ids = []
        list_values = []
        list_type_prescription = []
        list_dim = []
        list_centroids = []
        list_directions = []
        list_well_type = []
        list_model = []
        list_radius = []
        list_pbh = []
        list_flow_rate = []
        list_req = []
        list_wi = []
        list_block_dimension = []
        for id_text in cls.wells.keys():
            well = AllWells.wells[id_text]
            list_ids.append(well.id_text)
            list_volumes_ids.append(well.volumes_ids)
            list_values.append(well.value)
            list_type_prescription.append(well.type_prescription)
            list_dim.append(well.dim)
            list_centroids.append(well.centroids)
            list_directions.append(well.direction)
            list_well_type.append(well.well_type)
            list_model.append(well.model)
            list_radius.append(well.well_radius)
            list_pbh.append(well.pbh)
            list_flow_rate.append(well.flow_rate)
            list_req.append(well.req)
            list_wi.append(well.wi)
            list_block_dimension.append(well.block_dimension)

        df = pd.DataFrame({
            'ids': list_ids,
            'volumes_ids': list_volumes_ids,
            'values': list_values,
            'type_prescription': list_type_prescription,
            'dim': list_dim,
            'centroids': list_centroids,
            'direction': list_directions,
            'well_type': list_well_type,
            'model': list_model,
            'well_radius': list_radius,
            'pbh': list_pbh,
            'flow_rate': list_flow_rate,
            'req': list_req,
            'wi': list_wi,
            'block_dimension': list_block_dimension
        })

        df.to_csv(cls.database_name, index=False)

    @classmethod
    def load_wells_from_database(cls):
        df = pd.read_csv(cls.database_name)
        list_of_wells = []
        for index in df.index:
            line = df.loc[index]
            well = Well(id_text=str(line['ids']))
            volumes_ids = line['volumes_ids'].replace(' ', ',')
            volumes_ids2 = eval(volumes_ids)
            centroids = line['centroids'].replace('\n ', ',')
            centroids = centroids.replace('[ ', '[')
            centroids = centroids.replace('  ', ' ')
            centroids = centroids.replace(' ', ',')
            centroids2 = np.zeros((len(volumes_ids2), 3))
            exec('centroids2[:] = np.array(' + centroids + ')')
            value = line['values']
            type_prescription = line['type_prescription']
            dim = line['dim']
            direction = line['direction']
            well_type = line['well_type']
            model = line['model']
            radius = line=['well_radius']
            pbh = line['pbh']
            flow_rate = line['flow_rate']
            req = line['req']
            wi = line['wi']
            block_dimension = line['block_dimension']
            well.set_well(volumes_ids2, centroids2, block_dimension, type_prescription=type_prescription,
                          value=value, dim=dim, direction=direction, well_type=well_type,
                          model=model, well_radius=radius)
            well.pbh = pbh
            well.flow_rate = flow_rate
            well.req = req
            well.wi = wi


            list_of_wells.append(well)

        return list_of_wells


class Well(AllWells):
    type_prescription_list = ['pressure', 'flow_rate']
    type_direction_list = ['x', 'y', 'z']
    well_types = ['injector', 'producer']

    def __init__(self, id_text=''):
        test_id_text(id_text)
        if id_text == '':
            id_text = str(len(AllWells.wells))
        self.id_text = id_text

    def set_well(self, volumes_ids, centroids, block_dimension, type_prescription='pressure', value=0.0, dim=3, direction='x',
                 well_type='injector', model='peaceman', well_radius=0.5):
        '''

        :param volumes_ids:
        :param centroids:
        :param type_prescription:
        :param value:
        :param dim:
        :param direction:
        :param well_type:
        :param model:
            model in ['peaceman']
        :param well_radius: must be in the same unit of mesh
        :return:
        '''
        test_type_prescription(type_prescription)
        test_direction(direction)
        test_well_type(well_type)
        self.volumes_ids = volumes_ids.copy()
        self.value = value
        self.type_prescription = type_prescription
        self.dim = dim
        self.centroids = centroids.copy()
        self.direction = direction
        self.well_type = well_type
        self.model = model
        self._n_phases = 0
        self.well_radius = well_radius
        self.block_dimension = block_dimension.copy()
        if self.type_prescription == 'pressure':
            self.pbh = self.value
            self.flow_rate = None
        else:
            self.flow_rate = value
            self.pbh = None
        AllWells.wells[self.id_text] = weakref.proxy(self)

    def update_volumes(self, volumes_ids, centroids, block_dimension):
        self.volumes_ids = volumes_ids
        self.centroids = centroids
        self.block_dimension = block_dimension

    def update_n_phases(self, n_phases):
        self.n_phases = n_phases

    def update_mobilities(self, mobilities):
        self.mobilities = mobilities

    def update_req(self, req):
        self.req = req

    @property
    def centroids(self):
        return self._centroids

    @centroids.setter
    def centroids(self, value):
        n_cols = 3
        len_vols = len(self.volumes_ids)
        if value.shape[0] != len_vols:
            raise ValueError(f'lenght of centroids ({value.shape[0]}) different from volumes ({len_vols})')
        if value.shape[1] != n_cols:
            raise ValueError(f'n columns of centroids must be 3')
        self._centroids = value

    @property
    def n_phases(self):
        return self._n_phases

    @n_phases.setter
    def n_phases(self, value):
        if not isinstance(value, int):
            raise ValueError(f'n_phases must be integer type')
        self._n_phases = value

    @property
    def mobilities(self):
        return self._mobilities

    @mobilities.setter
    def mobilities(self, value):
        self.test_n_phases()
        ncols = value.shape[1]
        nrows = value.shape[0]
        if nrows != len(self.volumes_ids):
            raise ValueError('Number of columns must be equal to number of volumes in the well')
        if ncols != self.n_phases:
            raise ValueError('Number of rows must be equal to number of phases')
        self._mobilities = value

    @property
    def req(self):
        return self._req

    @req.setter
    def req(self, value):
        self._req = value

    @property
    def phases(self):
        return self._phases

    @phases.setter
    def phases(self, value):
        self.test_n_phases()
        if len(value) != self.n_phases:
            raise ValueError('Phases must be the same lenght of n_phases')
        for v in value:
            if not isinstance(v, str):
                raise NameError('Phase name must be string')
        self._phases = value

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, value):
        self.test_n_phases()
        self._rho = value

    @property
    def wi(self):
        return self._wi

    @wi.setter
    def wi(self, value):
        self._wi = value

    def calculate_req(self, **kwargs):
        '''

        :param kwargs:
            if self.model == 'peaceman':
                kwargs = [well_permeability]
        :return: equivalent r of well
        '''
        if self.model == 'peaceman':
            well_permeability = kwargs.get('well_permeability')
            if self.direction == 'z':
                k00 = well_permeability[:, 0, 0]
                k11 = well_permeability[:, 1, 1]
                h0 = self.centroids[:, 0]
                h1 = self.centroids[:, 1]
            else:
                # TODO: add more functionality
                raise NotImplementedError
            req = peaceman.get_req_vec(k00, k11, h0, h1)

        return req

    def calculate_WI(self, block_dimension, req, well_radius, well_perbeability, direction):
        if self.model == 'peaceman':
            pass
        else:
            # TODO: add more functionality
            raise NotImplementedError
        if direction == 'z':
            i1k = 0
            i2k = 1
            ih = 0
        else:
            # TODO: add more functionality
            raise NotImplementedError
        k00 = well_perbeability[:,i1k,i1k]
        k11 = well_perbeability[:,i2k,i2k]
        h = block_dimension[:, ih]

        wi = peaceman.get_WI_vec(h, k00, k11, req, well_radius)
        return wi

    def calculate_pbh(self):
        if self.type_prescription == 'pressure':
            pass
        else:
            # TODO: add more functionality
            raise NotImplementedError

    def calculate_flow_rate(self):
        if self.type_prescription == 'flow_rate':
            pass
        else:
            # TODO: add more functionality
            raise NotImplementedError

    def __del__(self):
        # print(f'\n Well {self.id_text} was deleted \n')
        AllWells.wells.pop(self.id_text)

    def __str__(self):
        return f'\nid: {self.id_text},' \
               f'\ntype prescription: {self.type_prescription},' \
               f'\nvalue: {self.value}\n'

    def test_n_phases(self):
        if self.n_phases == 0:
            raise ValueError('Update the number of phases')


def test_type_prescription(type_prescription):
    if type_prescription not in Well.type_prescription_list:
        raise NameError(f'type_prescription must be in {Well.type_prescription_list}')


def test_direction(direction):
    if direction not in Well.type_direction_list:
        raise NameError(f'direction must be in {Well.type_direction_list}')


def test_region(region):
    region_list = ['box']
    if region not in region_list:
        raise NameError(f'region must be in {region_list}')


def test_well_type(well_type):
    if well_type not in Well.well_types:
        raise NameError(f'well_type must be in {Well.well_types}')
