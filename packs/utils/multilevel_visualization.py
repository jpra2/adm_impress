from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh
import numpy as np

def get_faces_ids_by_level(mesh_level: np.ndarray, gids: np.ndarray, mesh_intersect_faces: np.ndarray, level):
    
        gid = gids[level]
        all_intersect_faces = mesh_intersect_faces[level]
        volumes = np.unique(gid[mesh_level == level])
        intersect_faces = np.unique(np.concatenate(all_intersect_faces[volumes]))
        return intersect_faces
        
    

def visualize_levels(M: FineScaleMesh, mesh_level: np.ndarray, gids: np.ndarray, mesh_faces: np.ndarray, name: str):
    
    ext = '.vtk'
    
    # meshset = m_object.core.mb.create_meshset()
    # m_object.core.mb.add_entities(m1, m_object.core.all_volumes)
    # m_object.core.mb.write_file(file_name, [m1])
    
    meshset = M.core.mb.create_meshset()
    faces_entities = M.core.all_faces
    volumes_entities = M.core.all_volumes
    print_function = M.core.print
    
    faces_to_export = np.ndarray([])
    
    all_levels = np.unique(mesh_level)
    
    for level in all_levels:
        faces_to_export = np.append(
            faces_to_export,
            get_faces_ids_by_level(mesh_level, gids, mesh_faces, level)
        )
    
    import pdb; pdb.set_trace()
    
    faces_to_export = np.unique(faces_to_export)
    M.core.mb.add_entities(meshset, M.core.all_volumes)
    M.core.mb.add_entities(meshset, faces_entities[faces_to_export])
    M.core.write_file(name + ext, [meshset])
    print('\n Success to export multilevel mesh \n')
    
        
    
    
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    