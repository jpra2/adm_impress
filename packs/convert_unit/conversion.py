from . import constants


class Conversion:

    def __init__(self, wells, data_impress):
        self.wells = wells
        self.data_impress = data_impress

    def convert_English_to_SI(self):

        '''
        converte as unidades do sistema americano para o SI
        '''
        k0 = constants.psi_to_Pa()
        k1 = constants.bbldia_to_m3seg()

        self.wells['values_p_ini'] *= k0
        self.wells['values_p'] *= k0
        self.wells['values_q'] *= k1

        k2 = constants.pe_to_m()

        self.data_impress['hs'] *= k2
        self.data_impress['volume'] *= k2**3
        self.data_impress['NODES'] *= k2
        self.data_impress['centroid_volumes'] *= k2
        self.data_impress['centroid_faces'] *= k2
        self.data_impress['centroid_edges'] *= k2
        self.data_impress['centroid_nodes'] *= k2

        k3 = constants.milidarcy_to_m2()

        self.data_impress['area'] *= k2**2
        self.data_impress['permeability'] *= k3
        self.data_impress['k_harm'] *= k3
        self.data_impress['dist_cent'] *= k2

        areas = self.data_impress['area']
        k_harm_faces = self.data_impress['k_harm']
        dist_cent = self.data_impress['dist_cent']

        self.data_impress['pretransmissibility'] = (areas*k_harm_faces)/dist_cent
