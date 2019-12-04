from . import constants

def convert_English_to_SI(M):

    '''
    converte as unidades do sistema americano para o SI
    '''
    k0 = constants.psi_to_Pa()
    k1 = constants.bbldia_to_m3seg()

    M.contours.datas['values_p_ini'] *= k0
    M.contours.datas['values_p'] *= k0
    M.contours.datas['values_q'] *= k1
    M.contours.update_values()

    k2 = constants.pe_to_m()

    M.data['hs'] *= k2
    M.data['volume'] *= k2**3
    M.data['NODES'] *= k2
    M.data['centroid_volumes'] *= k2
    M.data['centroid_faces'] *= k2
    M.data['centroid_edges'] *= k2
    M.data['centroid_nodes'] *= k2

    k3 = constants.milidarcy_to_m2()

    M.data['area'] *= k2**2
    M.data['permeability'] *= k3
    M.data['k_harm'] *= k3
    M.data['dist_cent'] *= k2

    areas = M.data['area']
    k_harm_faces = M.data['k_harm']
    dist_cent = M.data['dist_cent']

    M.data['pretransmissibility'] = (areas*k_harm_faces)/dist_cent

    M.data.update_variables_to_mesh()
    M.data.export_variables_to_npz()
