from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.examples.benchmarks_monophasic.cross.test_cross_1 import(
    create_primal_ids,
    export_primal_ids,
    create_dual_ids,
    export_dual_ids
)
import os
from packs.examples.benchmarks_monophasic.spe10.others import set_permeability
from packs.manager import MeshData

def get_properties():
    
    fine_mesh_path = os.path.join('mesh', 'spe10', 'layer_spe10_perturbada.msh')
    fine_mesh_properties_name = 'spe10_perturbada' 
    fine_mesh_path_v4 = fine_mesh_path

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)

    return fine_properties, fine_mesh_path


def run4():

    fp, fine_mesh_path = get_properties()
    set_permeability(fp, layer=1)
    
    bool_export_primal_id = True
    bool_export_dual_id = True
    my_dual_type = 1
    
    v1 = fp['permeability'][:,0,0]
    
    import pdb; pdb.set_trace()
    
    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('permx')
    mesh_data.insert_tag_data('permx', fp['permeability'][:,0,0], 'faces')
    mesh_data.export_all_elements_type_to_vtk('teste1_spe', 'faces')
    
    import pdb; pdb.set_trace()