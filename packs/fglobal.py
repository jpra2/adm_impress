import numpy as np
from packs.manager import SuperArrayManager


def local_names():
    my_dict = {
        'spe_data_name': 'spe_data'
    }
    
    return my_dict

def save_ijk_spe(spe_data: SuperArrayManager):
    volumes_centroids = spe_data['volumes_centroids']
    dxyz = np.array([[60.0, 220.0, 85.0]])/2
    
    divis = np.array([[60.0, 220.0, 85.0]])
    
    c2 = (volumes_centroids - dxyz)/divis
    
    ijk = c2.astype(np.int64)
    
    spe_data.insert_or_update_data({'ijk': ijk})
    spe_data.export_data()

def identify_spe_volume_by_ijk(spe_data: SuperArrayManager, local_ijk: np.ndarray):
    ijk_spe = spe_data['ijk']
    volumes_ids_spe = spe_data['volumes_ids']
    
    testi = ijk_spe[:, 0] = local_ijk[0]
    testj = ijk_spe[:, 1] = local_ijk[1]
    testk = ijk_spe[:, 2] = local_ijk[2]
    test = testi | testj | testk
    
    return volumes_ids_spe[test]

def identify_spe_volumes_by_ijks(spe_data: SuperArrayManager, ijks):
    volumes_ids = []
    for ijk in ijks:
        volumes_ids.append(
            identify_spe_volume_by_ijk(spe_data, ijk)
        )
    
    volumes_ids = np.concatenate(volumes_ids)
    return volumes_ids

