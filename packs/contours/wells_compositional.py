import numpy as np


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
        values_q = [] ## valor da vazao prescrita
        values_q_type = []
        q_composition = []
        ksi = []

        for p in data_wells:

            well = data_wells[p]
            type_region = well['type_region']
            tipo = well['type']
            prescription = well['prescription']
            value = (well['value']).astype(float)

            if type_region == direc.types_region_data_loaded[1]: #box

                p0 = well['p0']
                p1 = well['p1']
                limites = np.array([p0, p1])
                vols = get_box(centroids, limites)

                nv = len(vols)

                if prescription == 'Q':
                    val = value/nv
                    if tipo == 'Producer':
                        val *= -1
                    ws_q.append(vols)
                    values_type = np.repeat(well['value_type'], nv)
                    values = np.repeat(value, nv)
                    zs = np.repeat(well['z'], nv)
                    ksis = np.repeat(well['ksi_total'], nv)

                    ksi.append(ksis)
                    q_composition.append(zs)
                    values_q = np.append(values)
                    values_q_type.append(values_type)

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
        self['all_wells'] = np.union1d(ws_inj, ws_prod)
        self['values_p_ini'] = values_p.copy()
        self['value_type'] = values_q_type
        self['z'] = q_composition
        self['ksi_total'] = ksi
