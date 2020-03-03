from .. import directories as direc
import numpy as np



def set_saturation_regions(M, wells):

    centroids = M.data['centroid_volumes']
    n = len(centroids)

    for reg in direc.data_loaded[direc.names_data_loaded_lv0[4]]:
        d0 = direc.data_loaded[direc.names_data_loaded_lv0[4]][reg]
        tipo = d0['type']
        value = d0['value']

        if tipo == direc.types_region_for_saturation[0]:
            tamanho_variavel = len(M.data[M.data.variables_impress['saturation']])
            data = np.repeat(value, tamanho_variavel)
            M.data[M.data.variables_impress['saturation']] = data

        elif tipo == direc.types_region_for_saturation[1]:
            type1 = d0['type1_well']
            type2 = d0['type2_well']
            tipos = [type1, type2]
            all_wells = []
            all_values = []

            for tipo in tipos:
                if tipo == 'dirichlet':
                    wells = wells['ws_p']
                elif tipo == 'neumann':
                    wells = wells['ws_q']
                elif tipo == 'Injector':
                    wells = wells['ws_inj']
                elif tipo == 'Producer':
                    wells = wells['ws_prod']
                else:
                    continue

                all_wells.append(wells)
                all_values.append(np.repeat(value, len(wells)))
            

        else:
            raise("Tipo não suportado, disponíveis: all, box, wells")

            all_wells = np.array(all_wells)
            all_values = np.array(all_values)
            M.data[M.data.variables_impress['saturation']][all_wells] = all_values
