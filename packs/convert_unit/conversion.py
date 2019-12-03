from . import constants

def convert1(M):

    '''
    converte as unidades do sistema americano para o SI
    '''

    M.contours.datas['values_p_ini'] *= constants.psi_to_Pa()
    M.contours.datas['values_p'] *= constants.psi_to_Pa()
    M.contours.datas['values_q'] *= constants.bbldia_to_m3seg()
    M.data['hs'] *= constants.pe_to_m()
    M.data['volume'] *= constants.pe_to_m()**3
    M.data['NODES'] *= constants.pe_to_m()
    M.data['centroid_volumes'] *= constants.pe_to_m()
    M.data['centroid_faces'] *= constants.pe_to_m()
    M.data['centroid_edges'] *= constants.pe_to_m()
    M.data['centroid_nodes'] *= constants.pe_to_m()

    M.data['area'] *= constants.pe_to_m()**2
    M.data['permeability'] *= constants.milidarcy_to_m2()
    M.data['k_harm'] *= constants.milidarcy_to_m2()
    M.data['dist_cent'] *= constants.pe_to_m()

    areas = M.data['area']
    k_harm_faces = M.data['k_harm']
    dist_cent = M.data['dist_cent']

    M.data['pretransmissibility'] = (areas*k_harm_faces)/dist_cent

    M.data.update_variables_to_mesh()
    M.data.export_variables_to_npz()
