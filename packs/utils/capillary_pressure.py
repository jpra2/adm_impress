import numpy as np
from ..directories import data_loaded
import scipy.sparse as sp
from ..convert_unit import constants

class capillaryPressureBiphasic:
    names_pressure = ['psi', None, 'Pa']

    def __init__(self):
        '''
            Load the capillary_pressure.dat data in data file for interpolation
        '''
        self.cap_pressure = data_loaded['capillary_pressure']
        self.convert_capillary_pressure_to = data_loaded['convert_capillary_pressure_to']
        self.load_func()

    def load_func(self):
        if not self.cap_pressure:
            return 1

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
        return 0

    def get_pcow_from_sat(self, sat):
        '''
            get pcow from saturation on faces
        '''
        self.validate_arr(sat)
        if not self.cap_pressure:
            return np.zeros(len(sat))

        sat2 = np.array(sat)

        inds_menores = sat2 <= self.sat_min
        inds_maiores = sat2 >= self.sat_max
        inds_outros = ~(inds_maiores | inds_menores)

        resp = np.zeros(len(sat2))

        resp[inds_outros] = self.func(sat2[inds_outros])
        resp[inds_menores] = np.repeat(self.pcow_max, inds_menores.sum())
        resp[inds_maiores] = np.repeat(self.pcow_min, inds_maiores.sum())
        resp *= self.convert_cap_pressure()

        return resp

    def convert_cap_pressure(self):
        if self.convert_capillary_pressure_to == None:
            k = 1
        elif self.convert_capillary_pressure_to == 'psi':
            k = constants.atm_to_psi()
        elif self.convert_capillary_pressure_to == 'Pa':
            k = constants.atm_to_pa()
        else:
            raise NameError('Nome nao identificado')

        return k

    def validate_arr(self, val):

        if isinstance(val, list) or isinstance(val, tuple) or isinstance(val, np.ndarray):
            pass
        else:
            raise TypeError('sat deve ser um iteravel, menos um dict')

        if max(val) > 1.0 or min(val) < 0.0:
            raise ValueError(f'Saturacao errada. Min: {min(val)}. Max: {max(val)}')

    def get_gradient_cap_presure(self, dh, cap_pressure):
        return -(cap_pressure[:,1] - cap_pressure[:,0])/dh

    def get_flux_cap_faces(self, gradient_cap_faces, areas, k_harm, lambda_w_faces):
        return gradient_cap_faces*areas*k_harm

    def get_flux_cap_volumes(self, vols_viz_faces, flux_cap_faces):
        lines = np.concatenate([vols_viz_faces[:,0], vols_viz_faces[:,1]])
        n = len(np.unique(lines))
        cols = np.zeros(len(lines), dtype=int)
        data = np.concatenate([flux_cap_faces, -flux_cap_faces])
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n,1)).toarray().flatten()

        return flux_volumes

    def get_capillary_pressure_flux(self,
        faces: 'faces onde se vai calcular o gradiente da pressao capilar',
        vols_viz_faces: 'volumes vizinhos das faces',
        areas: 'areas das faces',
        k_harm: 'perm equivalente nas faces',
        dh: 'delta h',
        lambda_w: 'mobilidade de agua nas faces',
        cap_pressure: 'vetor com duas colunas: pcap dos vols_viz_faces'):

        if self.cap_pressure:
            # gradient_cap_faces = self.get_gradient_cap_presure(dh, cap_pressure)
            # flux_cap_faces = self.get_flux_cap_faces(gradient_cap_faces, areas, k_harm)
            flux_cap_faces = self.get_flux_cap_faces(self.get_gradient_cap_presure(dh, cap_pressure), areas, k_harm, lambda_w)
            flux_cap_volumes = self.get_flux_cap_volumes(vols_viz_faces, flux_cap_faces)
        else:
            flux_cap_volumes = np.zeros(len(np.unique(vols_viz_faces.flatten())))
            flux_cap_faces = np.zeros(len(faces))

        return flux_cap_faces, flux_cap_volumes
