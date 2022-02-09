import weakref
import os
from packs.directories import only_mesh_name
from packs.well_model import peaceman
import pandas as pd
import numpy as np
from scipy.sparse import find
import scipy.sparse as sp

from packs.properties.phisical_properties import PhisicalProperties
from packs.global_debug import GlobalDebug


class DebugPeacemanWellModel(GlobalDebug):
    # _debug = True
    pass


class AllWells:
    wells = dict()
    database_name = os.path.join('flying', 'wells_database_' + only_mesh_name + '.csv')
    # database_name_h5 = os.path.join('flying', 'wells_database_' + only_mesh_name + '.h5')

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
        list_z_well = []
        list_id_int = []
        list_n_phases = []
        list_rho_phases = []
        list_mobilities = []
        list_phases = []

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
            list_z_well.append(well.z_well)
            list_id_int.append(well.id_int)

            list_n_phases.append(well.n_phases)
            list_rho_phases.append(well.rho)
            list_mobilities.append(well.mobilities)
            list_phases.append(well.phases)

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
            'block_dimension': list_block_dimension,
            'z_well': list_z_well,
            'id_int' : list_id_int,
            'n_phases': list_n_phases,
            'rho_phases': list_rho_phases,
            'mobilities': list_mobilities,
            'phases': list_phases
        })

        df.to_csv(cls.database_name, index=False)

    @classmethod
    def load_wells_from_database(cls):
        df = pd.read_csv(cls.database_name)
        list_of_wells = []
        for index in df.index:
            line = df.loc[index]
            # well = Well(id_text=str(line['ids']))
            volumes_ids = line['volumes_ids'].replace(' ', ',')
            volumes_ids = volumes_ids.replace('[,', '[')
            volumes_ids = volumes_ids.replace(',,', ',')
            volumes_ids2 = np.array(eval(volumes_ids))
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
            radius = line['well_radius']
            pbh = line['pbh']
            flow_rate = line['flow_rate']
            req = line['req'].replace('\n', '')
            req = req.replace(' ', ',')
            req = np.array(eval(req))
            wi = line['wi'].replace('\n', '')
            wi = wi.replace(' ', ',')
            wi = np.array(eval(wi))
            block_dimension = line['block_dimension'].replace('\n ', ',')
            block_dimension = block_dimension.replace(' ', ',')
            block_dimension = np.array(eval(block_dimension))
            z_well = line['z_well']
            n_phases = line['n_phases']
            rho_phases = eval(line['rho_phases'])
            mobilities = line['mobilities'].replace('\n', '')
            mobilities = mobilities.replace(' ', ',')
            mobilities = np.array(eval(mobilities))
            phases = eval(line['phases'])

            well = Well()
            well.set_well(volumes_ids2, centroids2, block_dimension, type_prescription=type_prescription,
                          value=value, dim=dim, direction=direction, well_type=well_type,
                          model=model, well_radius=radius, id_text=str(line['ids']), z_well=z_well, n_phases=n_phases,
                          rho_phases=rho_phases, initialize=False)
            well.pbh = pbh
            well.flow_rate = flow_rate
            well.req = req
            well.wi = wi
            well.id_int = line['id_int']
            well.phases = phases
            well.mobilities = mobilities

            list_of_wells.append(well)


        return list_of_wells

    @classmethod
    def reorganize_id_int_wells(cls):
        for i, well in enumerate(AllWells.wells.values()):
            well.id_int = i

    @classmethod
    def update_wells_pbh(cls, wells, pressure1, volumes):
        """

        :param wells: list of all wells
        :param pressure1: pressure field with well pressures
        :param volumes: ids of all volumes of the mesh
        """
        i = 1


        for well in wells:
            if well.type_prescription == 'pressure':
                continue
            elif well.type_prescription == 'flow_rate':
                well_pressure = pressure1[volumes.max() + i]
                well.pbh = well_pressure
                i += 1
            else:
                raise NotImplementedError

    @classmethod
    def update_wells_flow_rate(cls, wells, flux_volumes):
        """

        :param wells: wells list
        :param flux_volumes: total flux volumes
        :param volumes: ids volumes
        """

        for well in wells:
            if well.type_prescription == 'flow_rate':
                continue
            elif well.type_prescription == 'pressure':
                well_flux = flux_volumes[well.volumes_ids].sum()
                well.flow_rate = well_flux
            else:
                raise NotImplementedError

    @classmethod
    def update_flux_phases(cls, wells, **kwargs):
        name = 'flux_phase_'
        w1: Well = wells[0]
        n_phases = w1.n_phases

        phases_flux = []

        for i in range(n_phases):
            phases_flux.append(kwargs.get(name + str(i)))

        phases_flux = np.array(phases_flux)
        total_flux = phases_flux.sum(axis=0)

        for well in wells:
            soma = well.mobilities.sum(axis=1)
            soma = soma.reshape(len(soma), 1)
            propo = well.mobilities / soma
            for i in range(n_phases):
                phases_flux[i][well.volumes_ids] = -total_flux[well.volumes_ids] * propo[:, i] + phases_flux[i][well.volumes_ids]

            # if well.well_type == 'producer':
            #     for i in range(n_phases):
            #         phases_flux[i][well.volumes_ids] = -total_flux[well.volumes_ids]*propo[:,i] + phases_flux[i][well.volumes_ids]
            # elif well.well_type == 'injector':
            #     for i in range(n_phases):
            #         phases_flux[i][well.volumes_ids] = -total_flux[well.volumes_ids]*propo[:,i] + phases_flux[i][well.volumes_ids]
            # else:
            #     raise NotImplementedError

        return phases_flux

    @classmethod
    def generate_report(cls, wells):
        pass













class Well:
    type_prescription_list = ['pressure', 'flow_rate']
    type_direction_list = ['x', 'y', 'z']
    well_types = ['injector', 'producer']

    def set_well(self, volumes_ids: np.ndarray, centroids: np.ndarray, block_dimension: np.ndarray, type_prescription='pressure', value=0.0, dim=3, direction='z',
                 well_type='injector', model='peaceman', well_radius=0.1, id_text='', z_well=1.0, n_phases=1,
                 well_permeability=None, rho_phases=None, initialize=True):
        '''

        :param volumes_ids:
        :param centroids:
        :param type_prescription: ['pressure', 'flow_rate']
        :param value:
        :param dim:
        :param direction: ['x', 'y', 'z']
        :param well_type: ['injector', 'producer']
        :param model: ['peaceman']
        :param well_radius: must be in the same unit of mesh dimension
        :return: None
        '''
        test_id_text(id_text)
        if id_text == '':
            id_text = str(len(AllWells.wells))
        self.id_text = id_text
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
        self._n_phases = n_phases
        self.mobilities = np.zeros((len(volumes_ids), n_phases))
        self.phases = ['phase_' + str(i) for i in range(n_phases)]
        self.well_radius = well_radius
        self.block_dimension = block_dimension.copy()
        self.z_well = z_well
        self.id_int = len(AllWells.wells)
        if self.type_prescription == 'pressure':
            self.pbh = self.value
            self.flow_rate = None
        else:
            if well_type == 'injector':
                self.flow_rate = -value
            elif well_type == 'producer':
                self.flow_rate = value
            else:
                raise NotImplementedError
            self.pbh = None

        if initialize is True:
            if well_permeability is None:
                print('\n Update the equivalent radius and Well Index \n')
                raise NotImplementedError
            else:
                self.req = self.calculate_req(well_permeability=well_permeability)
                self.wi = self.calculate_WI(self.block_dimension, self.req, self.well_radius, well_permeability, self.direction)

        if rho_phases is None:
            print('\n Update phase densities \n')
            self.rho = np.zeros(n_phases)
        else:
            self.rho = rho_phases

        AllWells.wells[self.id_text] = weakref.proxy(self)


    def update_volumes(self, volumes_ids, centroids):
        self.volumes_ids[:] = volumes_ids
        self.centroids[:] = centroids

    def update_volumes_ids(self, volumes_ids):
        self.volumes_ids[:] = volumes_ids

    def update_n_phases(self, n_phases):
        self.n_phases = n_phases

    def update_mobilities(self, mobilities):
        if self.well_type == 'injector':
            pass
        elif self.well_type == 'producer':
            self.mobilities = mobilities
        else:
            raise NotImplementedError

    def update_req(self, req):
        self.req = req

    def update_prescription(self, type_prescription='pressure', value=0.0):
        test_type_prescription(type_prescription)
        self.type_prescription = type_prescription
        self.value = value

    @property
    def centroids(self):
        return self._centroids

    @centroids.setter
    def centroids(self, value):
        n_cols = 3
        len_vols = len(self.volumes_ids)
        if value.shape[0] != len_vols:
            import pdb; pdb.set_trace()
            raise ValueError('lenght of centroids ({0}) different from volumes ({1})'.format((value.shape[0], len_vols)))
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

    @property
    def dZ(self):
        return self.z_well - self.centroids[:, 2]

    def calculate_req(self, well_permeability=None):
        '''

        :param kwargs:
            if self.model == 'peaceman':
                kwargs = [well_permeability: permeability of volumes in well]
        :return: equivalent radius of well
        '''
        if self.model == 'peaceman':
            if well_permeability is None:
                raise ValueError('well_permeability is None')
            # well_permeability = kwargs.get('well_permeability')
            if self.direction == 'z':
                k00 = well_permeability[:, 0, 0]
                k11 = well_permeability[:, 1, 1]
                h0 = self.block_dimension[:, 0]
                h1 = self.block_dimension[:, 1]
            else:
                # TODO: add more functionality
                raise NotImplementedError
        req = peaceman.get_req_vec(k00, k11, h0, h1)

        return DebugPeacemanWellModel.return_ones(req)

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

        return DebugPeacemanWellModel.return_ones(wi)

    def get_total_coefs(self):
        k = -1
        return k * self.wi * (self.mobilities.sum(axis=1))

    def get_source_gravity(self):
        g_acc = get_g_acc_z()
        source = np.zeros(len(self.volumes_ids))

        for phase in range(self.n_phases):
            source += self.wi * self.mobilities[:, phase] * (
                    self.rho[phase] * g_acc * (self.dZ))

        # if self.well_type == 'injector':
        #     fraction = self.fraction_mobility(self.mobilities)
        #     rhoMix = self.rho_mix(self.mobilities)
        # elif self.well_type == 'producer':
        #     fraction = self.fraction_mobility(vols_mobilities[self.volumes_ids])
        #     rhoMix = self.rho_mix(vols_mobilities[self.volumes_ids])
        #
        # source += self.wi * (rhoMix * g_acc * (self.dZ))
        # import pdb; pdb.set_trace()

        return source

    def get_coefs_mult_pbh(self):
        if self.type_prescription == 'pressure':
            coefs = self.get_total_coefs()
            return coefs * self.pbh
        else:
            raise ValidationError('this function is only for pressure restriction')

    @property
    def fraction_mobility(self):
        soma = self.mobilities.sum(axis=1)
        soma = soma.reshape(len(soma), 1)
        fraction = self.mobilities / soma
        return fraction

    @property
    def rho_mix(self):
        rhoMix = (self.wi * (self.fraction_mobility * self.rho).sum(axis=1)).sum() / self.wi.sum()
        # if self.well_type == 'injector':
        #     rhoMix = np.mean((self.fraction_mobility * self.rho).sum(axis=1))
        # elif self.well_type == 'producer':
        #
        #     rhoMix = (self.wi * self.fraction_mobility * self.rho).sum() / self.wi.sum()
        # else:
        #     raise NotImplementedError

        return rhoMix

    @property
    def pressure_drop(self):
        g_acc_z = get_g_acc_z()
        dp = g_acc_z * self.dZ * self.rho_mix
        return dp

    def __del__(self):
        # print(f'\n Well {self.id_text} was deleted \n')
        AllWells.wells.pop(self.id_text)
        AllWells.reorganize_id_int_wells()

    def __str__(self):
        return f'\nid: {self.id_text},' \
               f'\ntype prescription: {self.type_prescription},' \
               f'\nvalue: {self.value}\n'

    def test_n_phases(self):
        if self.n_phases == 0:
            raise ValueError('Update the number of phases')

    def report(self):
        pass







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


def test_id_text(id_text):
    if id_text in AllWells.wells.keys():
        raise NameError(f'id_text not must be in {AllWells.wells.keys()}')



def get_linear_problem(wells_list, T, gravity_source_term, volumes_ids=None):

    T2 = T.copy()
    q2 = gravity_source_term.copy()
    # lines = []
    # cols = []
    # data = []

    # import pdb; pdb.set_trace()
    for well in wells_list:
        if well.type_prescription == 'pressure':
            T2, q2 = insert_well_pressure_restriction(well, T2, q2)
            # T2, q2 = insert_well_pressure_restriction2(well, T2, q2)
            # id_wells.append(max_id_volumes+well.id_int)
            # pbh_list.append(well.pbh)

        elif well.type_prescription == 'flow_rate':
            T2, q2 = insert_well_flow_rate_prescription(well, T2, q2)


        else:
            # TODO: add more functionality
            raise NotImplementedError

    # inds_T = find(T2)
    #
    # lines = np.concatenate([inds_T[0], id_wells])
    # cols = np.concatenate([inds_T[1], id_wells])
    # data = np.concatenate([inds_T[2], np.ones(len(pbh_list))])
    # new_max_id = max_id_volumes + len(wells_list)
    # T2 = sp.csc_matrix((data, (lines, cols)), shape=(new_max_id, new_max_id))
    #
    # q3 = np.zeros(new_max_id)
    # q3[volumes_ids] = q2
    # q3[np.arange(volumes_ids.max()+1, new_max_id)] = pbh_list

    # return T2, q3

    return T2, q2


def get_g_acc_z():
    g_acc = PhisicalProperties._gravity_vector[2]
    return g_acc


def insert_well_pressure_restriction(well: Well, T: sp.csc_matrix, q: np.ndarray):
    """

    :param wells: well_object
    :param T: Transmissibility Matrix
    :param q: source term
    :param volumes_ids: all fine volumes ids
    :return:
    """
    coefs = well.get_total_coefs() # WI*sum(mobilities)
    T[well.volumes_ids, well.volumes_ids] += coefs
    source = well.get_coefs_mult_pbh() + well.get_source_gravity()
    q[well.volumes_ids] += source

    return T, q


def insert_well_pressure_restriction2(well: Well,T: sp.csc_matrix, q: np.ndarray):

    # coefs = -1 * well.wi * (well.mobilities.sum(axis=1))
    coefs = well.get_total_coefs() # = WI*sum(mobilities)
    T[well.volumes_ids, well.volumes_ids] += coefs
    # source = np.zeros(len(well.volumes_ids))
    #
    # for phase in range(well.n_phases):
    #     # source += 1 * well.wi * well.mobilities[:,phase] * (well.pbh + (well.rho[phase] * g_acc * (well.z_well - well.centroids[:,2])))
    #     source += 1 * well.wi * well.mobilities[:, phase] * (well.rho[phase] * g_acc * (well.z_well - well.centroids[:, 2]))

    source = well.get_source_gravity()

    q[well.volumes_ids] = source

    ff = find(T)
    n_gids = len(q) + 1
    gid_well = n_gids-1
    q2 = np.zeros(n_gids)
    T2 = sp.lil_matrix((n_gids, n_gids))

    T2[ff[0], ff[1]] = ff[2]
    repeated_gid_well = np.repeat(gid_well, len(well.volumes_ids))
    T2[well.volumes_ids, repeated_gid_well] = -coefs
    T2[gid_well, gid_well] = 1
    q2[:gid_well] = q
    q2[gid_well] = well.pbh

    return T2, q2


def insert_well_flow_rate_prescription(well: Well,T: sp.csc_matrix, q: np.ndarray):
    # g_acc = get_g_acc_z()

    # coefs = -1 * well.wi * (well.mobilities.sum(axis=1))
    # import pdb; pdb.set_trace()
    coefs = well.get_total_coefs()
    T[well.volumes_ids, well.volumes_ids] += coefs

    source = well.get_source_gravity()
    q[well.volumes_ids] += source

    ff = find(T)
    n_gids = len(q) + 1
    gid_well = n_gids - 1
    q2 = np.zeros(n_gids)
    T = sp.lil_matrix((n_gids, n_gids))

    T[ff[0], ff[1]] = ff[2]
    del ff
    repeated_gid_well = np.repeat(gid_well, len(well.volumes_ids))
    T[well.volumes_ids, repeated_gid_well] = -coefs
    T[repeated_gid_well, well.volumes_ids] = -coefs
    T[gid_well, gid_well] = coefs.sum()
    q2[:gid_well] = q
    q2[gid_well] = well.flow_rate - source.sum()

    return T.tocsc(), q2