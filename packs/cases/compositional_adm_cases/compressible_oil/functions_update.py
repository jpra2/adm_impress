import numpy as np

def update_global_vector_for_volumes_adjacencies_variable(global_vector_update, adjacencies, volume_variable, delta_lim):
    data = volume_variable[adjacencies]
    delta = np.abs(data[:,0] - data[:, 1])
    test = delta >= delta_lim
    if test.sum() > 0:  
        vols = np.unique(np.concatenate(adjacencies[test]))
        global_vector_update[vols] = True

def update_global_vector_for_latest_variable(global_vector_update, latest_variable, new_variable, delta_lim):
    delta = np.abs(latest_variable - new_variable)
    test = (delta >= delta_lim)
    latest_variable[test] = new_variable[test]
    global_vector_update[test.flatten()] = True