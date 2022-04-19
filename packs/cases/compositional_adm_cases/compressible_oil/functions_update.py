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
    # latest_variable[test] = new_variable[test]
    global_vector_update[test.flatten()] = True

def set_level0_delta_sat(level, saturation, adjacencies, delta_lim):
    delta_sat = np.abs(saturation[adjacencies[:, 1]] - saturation[adjacencies[:, 0]])
    test = delta_sat >= delta_lim
    if test.sum() > 0:
        vols = np.unique(np.concatenate(adjacencies[test]))
        level[vols] = 0

def update_global_vector_by_internal_face_variable(global_vector_update, latest_face_variable, new_face_variable, adjacencies, delta_lim):
    delta = np.abs(new_face_variable - latest_face_variable)
    test = delta >= delta_lim
    if test.sum() > 0:
        vols = np.unique(np.concatenate(adjacencies[test]))
        global_vector_update[vols] = True
        latest_face_variable[test] = new_face_variable[test]
