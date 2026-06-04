from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh

def run():
    fine_mesh_path = '32x32u_BL_D.msh'
    fine_mesh_name = '32BLD'
    
    fp = preprocess_mesh(mesh_name=fine_mesh_path, mesh_properties_name=fine_mesh_name)
    import pdb; pdb.set_trace()
    
    