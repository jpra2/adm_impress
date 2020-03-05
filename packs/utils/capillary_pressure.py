import numpy as np

class capillaryPressureBiphasic:

    def __init__(self):
        '''
            Load the capillary_pressure.dat data in data file for interpolation
        '''

        self.load_func()

    def load_func(self):

        points = np.loadtxt('data/capillary_pressure.dat')
        x = points[:,0]
        # self.x = x.copy()
        pcows = points[:,-1]
        # self.pcs = pcows.copy()
        grau = len(x)-1
        coefs = np.polyfit(x, pcows, grau)
        func = np.poly1d(coefs)

        self.pcow_max = pcows.max()
        self.pcow_min = pcows.min()
        self.sat_min = x[pcows == self.pcow_max][0]
        self.sat_max = x[pcows == self.pcow_min][0]
        self.func = func

    def get_pcow_from_sat(self, sat):
        '''
            get pcow from saturation
        '''

        self.validate_arr(sat)
        sat2 = np.array(sat)

        inds_menores = sat2 <= self.sat_min
        inds_maiores = sat2 >= self.sat_max
        inds_outros = ~(inds_maiores | inds_menores)

        resp = np.zeros(len(sat2))

        resp[inds_outros] = self.func(sat2[inds_outros])
        resp[inds_menores] = np.repeat(self.pcow_max, inds_menores.sum())
        resp[inds_maiores] = np.repeat(self.pcow_min, inds_maiores.sum())

        return resp

    def validate_arr(self, val):

        if isinstance(val, list) or isinstance(val, tuple) or isinstance(val, np.ndarray):
            pass
        else:
            raise TypeError('sat deve ser um iteravel, menos um dict')

        if max(val) > 1.0 or min(val) < 0.0:
            raise ValueError(f'Saturacao errada. Min: {min(val)}. Max: {max(val)}')
