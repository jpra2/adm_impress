import pdb
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh
import numpy as np

def get_faces_ids_by_level(mesh_level: np.ndarray, gids: np.ndarray, mesh_intersect_faces: np.ndarray, level):
    
    gid = gids[level]
    coarse_faces = mesh_intersect_faces[level]
    volumes = np.unique(gid[mesh_level == level])
    faces = np.unique(np.concatenate(coarse_faces[volumes]))
    return faces
    
    

def visualize_levels(M: FineScaleMesh, mesh_level: np.ndarray, gids: np.ndarray, mesh_faces: np.ndarray, name: str, loop: int, boundary_faces: np.ndarray):
    
    ext = '.vtk'
    
    # meshset = m_object.core.mb.create_meshset()
    # m_object.core.mb.add_entities(m1, m_object.core.all_volumes)
    # m_object.core.mb.write_file(file_name, [m1])
    M.data.update_variables_to_mesh()
    meshset = M.core.mb.create_meshset()
    meshset_faces = M.core.mb.create_meshset()
    faces_entities = M.core.all_faces
    
    faces_to_export = []
    
    all_levels = np.unique(mesh_level)
    
    for level in all_levels:
        faces_to_export.append(get_faces_ids_by_level(mesh_level, gids, mesh_faces, level))
    
    faces_to_export = np.unique(np.concatenate(faces_to_export))
    faces_to_export = np.setdiff1d(faces_to_export, boundary_faces)
    M.core.mb.add_entities(meshset, M.core.all_volumes)
    # M.core.mb.add_entities(meshset, faces_entities[faces_to_export])
    M.core.mb.add_entities(meshset_faces, faces_entities[faces_to_export])
    M.core.mb.write_file(name + '_volumes_'  + str(loop) + ext, [meshset])
    M.core.mb.write_file(name + '_faces_' + str(loop) + ext, [meshset_faces])
    
    print('\n Success to export multilevel mesh \n')
    # import pdb; pdb.set_trace()
    
        
    
    
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    