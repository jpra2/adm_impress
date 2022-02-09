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
        values_q_vol =[] ## valor da vazao prescrita em m3/s
        inj_cond = []
        z = []
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
                    z.append(well['z'])
                    z *= nv
                    inj_cond.append(well['injection_condition'])
                    inj_cond*=nv

                if prescription == 'Q':
                    val = value/nv * np.array(well['z'])
                    if tipo == 'Producer':
                        val *= -1

                    ws_q.append(vols)
                    value_type = well['value_type']
                    values = val

                    if value_type == 'volumetric':
                        if inj_cond[-1] == 'surface':
                            values = (1 - well['z'][-1]) * well['ksi_total'] * val
                            values[-1] = val[-1] * well['ksi_total']
                    values_q_vol.append(val.tolist())
                    values_q_vol *= nv


                    vals = np.repeat(values, nv)
                    if len(values_q)>0:
                        vals = (vals).reshape((len(well['z']),nv))
                        values_q = np.concatenate((values_q,vals), axis=1)
                        if inj_cond == 'surface':
                            values_q = (vals / well['ksi_total']).sum(axis=0)
                        values_q = np.concatenate((values_q, values_q), axis=0)

                    else:
                        values_q = np.append(values_q, vals).reshape((len(well['z']),nv))
                        if inj_cond == 'surface':
                            values_q = (values_q / well['ksi_total']).sum(axis=0)

                elif prescription == 'P':
                    val = value
                    ws_p.append(vols)
                    values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    ws_inj.append(vols)
                elif tipo == 'Producer':
                    ws_prod.append(vols)


        ws_q = np.array(ws_q).flatten()
        ws_p = np.array(ws_p).flatten()
        values_p = np.array(values_p).flatten()
        values_q = np.array(values_q)#.flatten()
        ws_inj = np.array(ws_inj).flatten()
        ws_prod = np.array(ws_prod).flatten()

        self['ws_p'] = ws_p.astype(int)
        self['ws_q'] = ws_q.astype(int)
        self['ws_inj'] = ws_inj.astype(int)
        self['ws_prod'] = ws_prod.astype(int)
        self['values_p'] = values_p
        self['values_q'] = values_q
        self['all_wells'] = np.union1d(ws_inj, ws_prod).astype(int)
        self['values_p_ini'] = values_p.copy()
        self['values_q_vol'] = np.array(values_q_vol).T
        self['inj_cond'] = np.array(inj_cond).flatten()
        self['z'] = np.array(z)
