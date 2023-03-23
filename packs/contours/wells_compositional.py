from .. import directories as direc
from ..utils.utils_old import get_box, getting_tag
from pymoab import types
import numpy as np
from ..data_class.data_manager import DataManager
import collections
from ..utils.utils_old import get_box
from .wells import Wells

class WellsCompositional(Wells):
    def get_wells(self):
        assert not self._loaded
        M = self.mesh

        data_wells = direc.data_loaded['Wells']
        centroids = M.data['centroid_volumes']
        gravity = direc.data_loaded['gravity']


        ws_p = [] ## pocos com pressao prescrita
        ws_q = [] ## pocos com vazao prescrita
        ws_inj = [] ## pocos injetores
        ws_prod = [] ## pocos produtores
        values_p = [] ## valor da pressao prescrita
        values_q = np.array([]) ## valor da vazao prescrita
        values_qa = np.array([]) ## valor da vazao prescrita
        values_qco2 = np.array([]) ## valor da vazao prescrita
        values_q_vol =[] ## valor da vazao prescrita em m3/s
        values_q_vol_a =[]
        values_q_vol_co2 =[]
        inj_cond = []
        z = []
        za = []
        zco2 = []
        DT = []

        for p in data_wells:

            well = data_wells[p]
            type_region = well['type_region']
            tipo = well['type']
            prescription = well['prescription']
            value = np.array(well['value']).astype(float)

            if type_region == direc.types_region_data_loaded[1]: #box
                p0 = well['p0']
                p1 = well['p1']
                limites = np.array([p0, p1])
                vols = get_box(centroids, limites)

                nv = len(vols)
                if tipo == 'Injector':
                    value_a = np.array(well['value_a']).astype(float)
                    value_co2 = np.array(well['value_co2']).astype(float)
                    t_max = well['t_max']
                    nP_inj = well['nP_inj']
                    t_cicle = t_max/nP_inj
                    DT = np.linspace(t_cicle, t_max, num = nP_inj)
                    DT *= nv
                    z.append(well['z'])
                    z *= nv
                    za.append(well['za'])
                    za *= nv
                    zco2.append(well['zco2'])
                    zco2 *= nv
                    inj_cond.append(well['injection_condition'])
                    inj_cond*=nv
                    #import pdb; pdb.set_trace()

                if prescription == 'Q':
                    val = value/nv * np.array(well['z'])
                    val_a = value_a/nv * np.array(well['za'])
                    val_co2 = value_co2/nv * np.array(well['zco2'])


                    if tipo == 'Producer':
                        val *= -1

                    ws_q.append(vols)
                    value_type = well['value_type']
                    values = val
                    values_a = val_a
                    values_co2 = val_co2


                    '''if value_type == 'volumetric':
                        if inj_cond[-1] == 'surface':
                            values = (1 - well['z'][-1]) * well['ksi_total'] * val
                            values[-1] = val[-1] * well['ksi_total']'''

                    values_q_vol.append(val.tolist())
                    values_q_vol *= nv

                    values_q_vol_a.append(val_a.tolist())
                    values_q_vol_a *= nv

                    values_q_vol_co2.append(val_co2.tolist())
                    values_q_vol_co2 *= nv

                    vals = np.repeat(values, nv)
                    vals_a = np.repeat(values_a, nv)
                    vals_co2 = np.repeat(values_co2, nv)

                    if len(values_q)>0:
                        vals = (vals).reshape((len(well['z']),nv))
                        values_q = np.concatenate((values_q,vals), axis=1)
                        if inj_cond == 'surface':
                            values_q = (vals).sum(axis=0)
                            values_q = np.concatenate((values_q, values_q), axis=1)

                    else:
                        values_q = np.append(values_q, vals).reshape((len(well['z']),nv))
                        if inj_cond == 'surface':
                            values_q = (values_q).sum(axis=0)

                    if len(values_qa)>0:
                        vals_a = (vals_a).reshape((len(well['za']),nv))
                        values_qa = np.concatenate((values_qa,vals_a), axis=1)
                        if inj_cond == 'surface':
                            values_qa = (vals_a).sum(axis=0)
                            values_qa = np.concatenate((values_qa, values_qa), axis=1)

                    else:
                        values_qa = np.append(values_qa, vals_a).reshape((len(well['za']),nv))
                        if inj_cond == 'surface':
                            values_qa = (values_qa).sum(axis=0)

                    if len(values_qco2)>0:
                        vals_co2 = (vals_co2).reshape((len(well['zco2']),nv))
                        values_qco2 = np.concatenate((values_qco2,vals_co2), axis=1)
                        if inj_cond == 'surface':
                            values_qco2 = (vals_co2).sum(axis=0)
                            values_qco2 = np.concatenate((values_qco2, values_qco2), axis=1)

                    else:
                        values_qco2 = np.append(values_qco2, vals_co2).reshape((len(well['zco2']),nv))
                        if inj_cond == 'surface':
                            values_qco2 = (values_qco2).sum(axis=0)

                elif prescription == 'P':
                    val = value
                    ws_p.append(vols)
                    values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    ws_inj.append(vols)
                elif tipo == 'Producer':
                    ws_prod.append(vols)
                #import pdb; pdb.set_trace()
        ws_q = np.array(ws_q).flatten()
        ws_p = np.array(ws_p).flatten()
        values_p = np.array(values_p).flatten()
        values_q = np.array(values_q)#.flatten()
        values_qa = np.array(values_qa)#.flatten()
        values_qco2 = np.array(values_qco2)#.flatten()
        ws_inj = np.array(ws_inj).flatten()
        #import pdb; pdb.set_trace()
        ws_prod = np.array(ws_prod).flatten()
        self['ws_p'] = ws_p.astype(int)
        self['ws_q'] = ws_q.astype(int)
        self['ws_inj'] = ws_inj.astype(int)
        self['ws_prod'] = ws_prod.astype(int)
        self['values_p'] = values_p
        self['values_q'] = values_q
        self['values_qa'] = values_qa
        self['values_qco2'] = values_qco2
        self['all_wells'] = np.union1d(ws_inj, ws_prod).astype(int)
        self['values_p_ini'] = values_p.copy()
        self['values_q_vol'] = np.array(values_q_vol).T
        self['values_q_vol_a'] = np.array(values_q_vol_a).T
        self['values_q_vol_co2'] = np.array(values_q_vol_co2).T
        self['inj_cond'] = np.array(inj_cond).flatten()
        self['DT'] = np.array(DT)
        self['z'] = np.array(z)
        self['za'] = np.array(za)
        self['zco2'] = np.array(zco2)
