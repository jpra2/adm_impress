from packs.manager.meshmanager import MeshProperty
import numpy as np
import os

def set_permeability(fp:MeshProperty, layer=1):
    
    perm_path = os.path.join('data', f'layer_{layer}', 'perm_layer.npy')
    perm = np.load(perm_path)
    perm = perm[:, 0:2, 0:2]
    
    fp.insert_or_update_data({'permeability': perm})
    fp.export_data()
    
    